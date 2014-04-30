/*
 * LocalLuscher.hpp
 *
 *  Created on: Jun 20, 2013
 *      Author: Thibaut Metivet
 */

#ifndef LOCAL_LUSCHER_HPP_
#define LOCAL_LUSCHER_HPP_

#include "LQCDA.hpp"

namespace LocalLuscher {

    typedef LQCDA::ResampledFitData<LQCDA::BootstrapResampler, std::vector<double>, std::vector<double>, LQCDA::FullCovariance, LQCDA::FullCovariance, LQCDA::FullCovariance> RSFitData;
    typedef LQCDA::FitDataWrapper<std::vector<double>, std::vector<double>> RSDataWrapper;
    //typedef LQCDA::ModelFunction<RSFitData::data_type, RSFitData::x_type> ModelFcn;
    typedef LQCDA::FitModel<std::vector<double>, std::vector<double>> Model;

} // namespace LocalLuscher

#endif	// LOCAL_LUSCHER_HPP_
