
#ifndef GUARD_dummyutilities_h
#define GUARD_dummyutilities_h
#include "tree.h"
#include "treefuns.h"
#include "bartfuns.h"
#include "info.h"
#include "Mvn.h"
#include <limits>
#include <iostream>
#include <fstream>
#include <map>
#include <tuple>
#include <algorithm>
#include <igl/slice.h>

namespace dummyutilities {
	//--------------------------------------------------
	typedef std::vector<int> index_v;
	typedef std::vector<index_v> index_2d_v;

	template <class T> void printVec(std::vector<T> const &input);
	template <class T> void printVec(std::vector<std::vector<T> > const &input);

	void VectorToEvector(index_v& v, Eigen::VectorXi& ev);
	void VectorToEvector(std::vector<double>& v, Eigen::VectorXd& ev);
	void EvectorToVector(Eigen::VectorXi& ev, index_v& v);
	void EvectorToVector(Eigen::VectorXd& ev, std::vector<double>& v);
	void EmatrixToVectorDouble(Eigen::MatrixXd& emat ,Eigen::VectorXi& ev);
	void parse2DCsvFile(std::string inputFileName, std::vector<std::vector<double>>& data);
	void readCSV(std::string file, int rows, int cols, Eigen::MatrixXd& data);
	int csvWrite(const Eigen::MatrixXd& inputMatrix, const std::string& fileName, const std::streamsize dPrec);
	int csvWrite(double* data, int n, const std::string& fileName, const std::streamsize dPrec);
	void moveItemToBack(index_2d_v& v2d, size_t itemIndex);
	void moveItemToBack(index_v& v1d, size_t itemIndex);
	void EMatToMap(Eigen::MatrixXd& EM, bool birth, size_t bdId, index_v& ordered_nids, std::map<std::tuple<size_t, size_t>, double>& mapA);
	//--------------------------------------------------
	// Tree related
	//--------------------------------------------------
	//Split the data ids for the bottom node
	void splitBnodeXids(index_v& xids, size_t v, size_t c, xinfo& xi, dinfo& di, index_2d_v& new_split_order);
	//--------------------------------------------------
	// Overload (treefuns.h): fit function for reordering
	void fit(tree& t, xinfo& xi, size_t p, size_t n, double *x,  double* fv, index_v& xid_nid);
	//--------------------------------------------------
	//Get ordered (ascending by nid) nid vector
	void getOrderedNid(tree::npv& bns, index_v& ordered_nids);
	//--------------------------------------------------
	//Reordering the data by tree fitting
	void reorder(index_v& xid_nid, index_v& ordered_nids, index_2d_v& new_order);
	//
	//--------------------------------------------------
	// Eigen related
	//--------------------------------------------------
	// Construct A from eA (the map of A entries)
	void getMatAfromeA(std::map<std::tuple<size_t, size_t>, double>& eA, index_v& ordered_nids, Eigen::MatrixXd& A);
	//Get matrix A
	void getMatA(pinfo& pi, index_2d_v& new_order, const Eigen::MatrixXd& SigmaInv, Eigen::MatrixXd& A);
	//Get matrix A block
	void getMatA(index_2d_v& new_order, const Eigen::MatrixXd& SigmaInv, Eigen::MatrixXd& Ab, index_2d_v& block_order);
	//--------------------------------------------------
	//Get scalar u
	double getScalarU(index_2d_v& new_order, const Eigen::MatrixXd& SigmaInv, const Eigen::MatrixXd& B,  double *r);
	//--------------------------------------------------
	//Get birth marginal likelihood ratio given T
	double getLogMglrBirth(tree& x, pinfo& pi, double *r, const Eigen::MatrixXd& SigmaInv, index_2d_v& new_order, index_2d_v& new_split_order, size_t birth_node_order, index_v& ordered_nids, bool sigmaNoChange, Eigen::MatrixXd& mA_i1);
	//--------------------------------------------------
	//Get death marginal likelihood ratio given T
	double getLogMglrDeath(tree& x, pinfo& pi, double *r, const Eigen::MatrixXd& SigmaInv, index_2d_v& new_order, size_t death_nodes_order, index_v& ordered_nids, bool sigmaNoChange, Eigen::MatrixXd& mA_i1);

	//--------------------------------------------------
	// Tree & Eigen related
	//--------------------------------------------------
	// Overload (bd.h): bd function for Metropolis-Hasting MCMC
	bool bd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double *r, const Eigen::MatrixXd& SigmaInv,
		std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen, index_v& xid_nid, bool sigmaNoChange);
	//--------------------------------------------------
	// Overload (bartfuns.h) : drmu for drawing dependent mus
	void drmu(tree& t, pinfo& pi, double *r, const Eigen::MatrixXd& SigmaInv, index_v& xid_nid, bool sigmaNoChange);

}

#endif