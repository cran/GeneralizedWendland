// [[Rcpp::depends(Rcpp, RcppEigen, BH)]]
#include "Wendland.h"

WendlandParam::WendlandParam(){
  checkParameters();
}

WendlandParam::WendlandParam(double _rng, double _psil, double _kap, double _mu,
  double _nug) : range(_rng), sill(_psil), kappa(_kap), mu(_mu), nugget(_nug){
  checkParameters();
}

void WendlandParam::checkParameters(){
  if (range < 0.0) Rcpp::stop("Range must not be negative.");
  if (sill < 0.0) Rcpp::stop("Sill must not be negative.");
  if (kappa < 0.0) Rcpp::stop("Kappa must not be negative.");
  if (mu < 0.0) Rcpp::stop("Mu must not be negative.");
  if (nugget < 0.0) Rcpp::stop("Nugget must not be negative.");

  if (sill + nugget == 0.0) Rcpp::stop("Produces zero valued covariance matrix.");
  if (mu < (1.0 + 2.0)/2.0 + kappa) Rcpp::warning("Mu < lambda(d, kappa). Covariance matrix may not be pd.");
}

Wendland::Wendland(){
  computeBetaConstant();
};

Wendland::~Wendland(){
  deleteInterpolator();
  deleteIntegrator();
}

void Wendland::computeBetaConstant(){
  beta_constant = boost::math::beta(1.0 + 2.0 * param.kappa, param.mu);
}

double Wendland::computeIntegral(double distance){
  double r = distance/param.range;
  return integrator->integrate(r, 1.0, [&](double x){
    return std::pow(x*x - r*r, param.kappa) * std::pow(1.0 - x, param.mu - 1.0);
  });
}

void Wendland::setParameters(double range, double sill, double kappa, double mu, double nugget){
  param = WendlandParam(range, sill, kappa, mu, nugget);
  computeBetaConstant();
  if (interpolator) interpolator->initialize();
}

void Wendland::setEpsTol(double _epstol){
  epstol = _epstol;
}

void Wendland::setIntegrator(double abstol, double reltol, int intervals=0, int qag_key=0){
  deleteIntegrator();
  integrator = new Integrator(abstol, reltol, intervals, qag_key);
}

void Wendland::deleteIntegrator(){
  if (integrator) {
    delete integrator;
    integrator = nullptr;
  }
}

void Wendland::setInterpolator(int num_points, int interp_type=0){
  deleteInterpolator();
  if (interp_type){
    interpolator = new Interpolator<Wendland>(this, interp_type, num_points, param.range);
    interpolator->initialize();
  }
}

void Wendland::deleteInterpolator(){
  if (interpolator) {
    delete interpolator;
    interpolator = nullptr;
  }
}

double Wendland::compute(const double& distance){
  if (!interpolator || !interpolator->isInitialized()){
    return (distance < param.range) ? computeIntegral(distance) : 0.0;
  } else {
    return (distance < param.range) ? interpolator->compute(distance) : 0.0;
  }
}


Rcpp::NumericVector Wendland::computeVector(const Rcpp::NumericVector& dvec){
  double d;
  int len = dvec.length();
  Rcpp::NumericVector cvec(len);

  for (int i = 0; i < len; i++){
    d = dvec(i);

    if (d>=param.range){
      continue;

    } else if (d<epstol){
      cvec(i)=param.sill+param.nugget;

    } else {
      cvec(i)=param.sill*compute(d)/beta_constant;
    }
  }
  return cvec;
}


Eigen::MatrixXd Wendland::computeMatrix(const Eigen::MatrixXd& dmat){
  double d;
  int nrow = dmat.rows();
  int ncol = dmat.cols();
  bool is_square = (nrow == ncol);
  Eigen::MatrixXd cmat(nrow, ncol);
  cmat.setZero();

  for (int i = 0; i < nrow; i++){
    for (int j = 0; j < ncol; j++){
      d = dmat(i,j);

      if (d >= param.range || cmat(i,j)){
        continue;

      } else if (is_square && cmat(j,i) && (d==dmat(j,i))){
        cmat(i,j) = cmat(j,i);

      } else if (d < epstol){
        cmat(i,j) = param.sill+(i==j?param.nugget:0.0);

      } else {
        cmat(i,j) = param.sill*compute(d)/beta_constant;
      }
    }
  }

  return cmat;
}


Eigen::SparseMatrix<double> Wendland::computeMSparse(const Eigen::SparseMatrix<double>& dmat){
  double d, value;
  int row, col;
  int nrow = dmat.rows();
  int ncol = dmat.cols();
  bool is_square = (nrow == ncol);

  Eigen::SparseMatrix<double> cmat(nrow, ncol);
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> computed(nrow, ncol);
  computed.setZero();

  std::vector<Eigen::Triplet<double>> tripletList;
  tripletList.reserve(nrow * ncol);

  for (int k=0; k<dmat.outerSize(); ++k){
    for (Eigen::SparseMatrix<double>::InnerIterator it(dmat,k); it; ++it){
      d = it.value();
      row = it.row();
      col = it.col();

      if ((d<param.range) && (!computed(row,col))){
        value = (d<epstol && row == col) ? param.sill+param.nugget :
        param.sill*compute(d)/beta_constant;

        if (value > epstol){
          tripletList.push_back(Eigen::Triplet<double>(row,col,value));
          computed(row,col)=1;

          if (is_square && !computed(col,row) && (d == dmat.coeff(col,row))){
            tripletList.push_back(Eigen::Triplet<double>(col,row,value));
            computed(col,row)=1;
          }
        }
      } else {
        continue;
      }
    }
  }

  tripletList.shrink_to_fit();
  cmat.setFromTriplets(tripletList.begin(), tripletList.end());
  return cmat;
}

Rcpp::List Wendland::computeSpam(const Eigen::MatrixXi& index, const Eigen::VectorXd& values){
  double value;
  int counter = 0;
  std::tuple<int,int,double> key, tkey;
  std::map<std::tuple<int,int,double>, double> map;

  for (int r=0; r<index.rows();r++){
    key = std::tuple<int,int,double>(index(r,0),index(r,1),values(r));
    tkey = std::tuple<int,int,double>(index(r,1),index(r,0),values(r));

    if ((std::get<2>(key) >= param.range) || (map.find(key) != map.end())){
      continue;

    } else {

      if (map.find(tkey) != map.end()){
        value = map.find(tkey)->second;

      } else {
        value = (std::get<2>(key) < epstol &&
          std::get<0>(key) == std::get<1>(key)) ? param.sill+param.nugget :
          param.sill*compute(std::get<2>(key))/beta_constant;
      }

      if (value > epstol){
        map.insert(std::pair<std::tuple<int,int,double>, double>(key, value));
      }
    }
  }

  Rcpp::IntegerMatrix new_index(map.size(),2);
  Rcpp::NumericVector cov_values(map.size());


  for (std::map<std::tuple<int,int,double>,double>::iterator it = map.begin(); it != map.end(); it++){
    Rcpp::IntegerMatrix::Row ind_ref = new_index(counter, Rcpp::_ );
    ind_ref = Rcpp::IntegerVector::create(std::get<0>(it->first), std::get<1>(it->first));
    cov_values(counter)= it->second;
    counter++;
  }

  Rcpp::List result = Rcpp::List::create(Rcpp::Named("indices") = new_index,
                                         Rcpp::Named("values") = cov_values);
  return result;
}

//void wendlandFinalizer(Wendland* object){
//  delete object;
//}

using namespace Rcpp;
RCPP_EXPOSED_CLASS(Wendland)
RCPP_MODULE(Wendland) {
  class_<Wendland>("Wendland")
  .constructor("Default QNG integration")
  .method("setParameters", &Wendland::setParameters)
  .method("setEpsTol", &Wendland::setEpsTol)
  .method("setIntegrator", &Wendland::setIntegrator)
  .method("setInterpolator", &Wendland::setInterpolator)
  .method("compute", &Wendland::compute)
  .method("computeVector", &Wendland::computeVector)
  .method("computeMatrix", &Wendland::computeMatrix)
  .method("computeSparse", &Wendland::computeMSparse)
  .method("computeSpam", &Wendland::computeSpam)
//  .finalizer(&wendlandFinalizer)
;}
