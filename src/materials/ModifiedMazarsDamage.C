/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*                       BlackBear                              */
/*                                                              */
/*           (c) 2017 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "ModifiedMazarsDamage.h"

#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "ElasticityTensorTools.h"
#include "MooseUtils.h"

#include "libmesh/utility.h"

registerMooseObject("BlackBearApp", ModifiedMazarsDamage);

InputParameters
ModifiedMazarsDamage::validParams()
{
  InputParameters params = ScalarDamageBase::validParams();
  params.addClassDescription("Mazars scalar damage model");
  params.addRequiredCoupledVar("tensile_strength",
                               "Tensile stress threshold for damage initiation");
  params.addRequiredCoupledVar("compressive_yield_strength",
                               "Compressive stress threshold for damage initiation");                               
  params.addRequiredParam<Real>("a_t",
                                "A_t parameter that controls the shape of the response in tension");
  params.addRequiredParam<Real>("b_t",
                                "B_t parameter that controls the shape of the response in tension");
  params.addRequiredParam<Real>(
      "a_c", "A_c parameter that controls the shape of the response in compression");
  params.addRequiredParam<Real>(
      "b_c", "B_c parameter that controls the shape of the response in compression");
  return params;
}

ModifiedMazarsDamage::ModifiedMazarsDamage(const InputParameters & parameters)
  : ScalarDamageBase(parameters),
    GuaranteeConsumer(this),
    _tensile_strength(coupledValue("tensile_strength")),
    _compressive_yield_strength(coupledValue("compressive_yield_strength")),    
    _a_t(getParam<Real>("a_t")),
    _a_c(getParam<Real>("a_c")),
    _b_t(getParam<Real>("b_t")),
    _b_c(getParam<Real>("b_c")),
    _kappa(declareProperty<Real>(_base_name + "kappa")),
    _kappa_old(getMaterialPropertyOld<Real>(_base_name + "kappa")),
    _stress(getMaterialProperty<RankTwoTensor>(_base_name + "stress")),
    _mechanical_strain(getMaterialProperty<RankTwoTensor>(_base_name + "mechanical_strain")),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),
    _eigval(3, 0.0),
    _positive_stress(3, 0.0),
    _positive_strain(3, 0.0)
{
}

void
ModifiedMazarsDamage::initQpStatefulProperties()
{
  ScalarDamageBase::initQpStatefulProperties();
  _kappa[_qp] = 0.0;
}

void
ModifiedMazarsDamage::initialSetup()
{
  if (!hasGuaranteedMaterialProperty(_elasticity_tensor_name, Guarantee::ISOTROPIC))
    mooseError("MazarsDamage requires that the elasticity tensor be guaranteed isotropic");
}

void
ModifiedMazarsDamage::updateQpDamageIndex()
{
  _stress[_qp].symmetricEigenvalues(_eigval);
  for (unsigned int i = 0; i < 3; ++i)
    _positive_stress[i] = std::max(_eigval[i], 0.0);

  _mechanical_strain[_qp].symmetricEigenvalues(_eigval);
  for (unsigned int i = 0; i < 3; ++i)
    _positive_strain[i] = std::max(_eigval[i], 0.0);

  const Real abs_eig_strain_1 = std::abs(_eigval[0]);
  const Real abs_eig_strain_3 = std::abs(_eigval[2]);
  const Real maxAbsValue_eig_strain = std::max(abs_eig_strain_1, abs_eig_strain_3);

  Real sign_max_strain;

  if (maxAbsValue_eig_strain == abs_eig_strain_1) {
      sign_max_strain = (_eigval[0] >= 0) ? 1 : -1;
  } else {
      sign_max_strain = (_eigval[3] >= 0) ? 1 : -1;
  }

  const Real equiv_tensile_strain =
      std::sqrt(Utility::pow<2>(_positive_strain[0]) + Utility::pow<2>(_positive_strain[1]) +
                Utility::pow<2>(_positive_strain[2]));

  Real alpha_t = 0.0;
  Real alpha_c = 0.0;

  const Real e = ElasticityTensorTools::getIsotropicYoungsModulus(_elasticity_tensor[_qp]);
  const Real nu = ElasticityTensorTools::getIsotropicPoissonsRatio(_elasticity_tensor[_qp]);

  const Real sum_pos_stress = _positive_stress[0] + _positive_stress[1] + _positive_stress[2];

  for (unsigned int i = 0; i < 3; ++i)
  {
    const Real epsilon_ti = ((1.0 + nu) * _positive_stress[i] - nu * sum_pos_stress) / e;
    const Real epsilon_ci = _eigval[i] - epsilon_ti; //_eigval contains the principal strains
    if (!MooseUtils::absoluteFuzzyEqual(equiv_tensile_strain, 0.0))
    {
      const Real ets2 = Utility::pow<2>(equiv_tensile_strain);
      alpha_t += epsilon_ti * _positive_strain[i] / ets2;
      alpha_c += epsilon_ci * _positive_strain[i] / ets2;
    }
  }

  const Real kappa_0t = _tensile_strength[_qp] / e; // Threshold strain for damage
  const Real kappa_0c = _compressive_yield_strength[_qp] / e; // Threshold strain for damage

  Real & kappa = _kappa[_qp];

  if (sign_max_strain>0)
    kappa = std::max(kappa_0t, _kappa_old[_qp]);
  else
    kappa = std::max(kappa_0c, _kappa_old[_qp]);

  _damage_index[_qp] = 0.0;
  if (equiv_tensile_strain > kappa)
  {
    kappa = equiv_tensile_strain;
    const Real tensile_damage =
        1.0 - kappa_0t * (1.0 - _a_t) / kappa - _a_t * std::exp(-_b_t * (kappa - kappa_0t));
    const Real compressive_damage =
        1.0 - kappa_0c * (1.0 - _a_c) / kappa - _a_c * std::exp(-_b_c * (kappa - kappa_0c));
    _damage_index[_qp] = std::max(alpha_t * tensile_damage + alpha_c * compressive_damage, 0.0);
  }

  // Prevevscode-file://vscode-app/Applications/Visual%20Studio%20Code.app/Contents/Resources/app/out/vs/code/electron-sandbox/workbench/workbench.htmlnt damage index from decreasing or from being greater than 1
  _damage_index[_qp] = std::max(_damage_index[_qp], _damage_index_old[_qp]);
  _damage_index[_qp] = std::min(_damage_index[_qp], 1.0);
}
