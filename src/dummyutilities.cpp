
#include "dummyutilities.h"

namespace dummyutilities{

	template <class T>
	void printVec(std::vector<T> const &input)
	{
	  for (int i = 0; i < input.size(); i++)
	  {
	    std::cout << input.at(i) << " ";
	  }

	  std::cout << std::endl;
	}


	template <class T>
	void printVec(std::vector<std::vector<T> > const &input)
	{
	  int n = input.size();
	  for (int i = 0; i < n; i++)
	  {
	    printVec(input[i]);
	  }

	  std::cout << std::endl;
	}

	//--------------------------------------------------
	void VectorToEvector(index_v& v, Eigen::VectorXi& ev)
	{
		ev = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(v.data(), v.size());
	}
	//--------------------------------------------------
	void VectorToEvector(std::vector<double>& v, Eigen::VectorXd& ev)
	{
		ev = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(v.data(), v.size());
	}
	//--------------------------------------------------
	void EvectorToVector(Eigen::VectorXi& ev, index_v& v)
	{
		v.resize(ev.size());
		Eigen::VectorXi::Map(&v[0], ev.size()) = ev;
	}
	//--------------------------------------------------
	void EvectorToVector(Eigen::VectorXd& ev, std::vector<double>& v)
	{
		v.resize(ev.size());
		Eigen::VectorXd::Map(&v[0], ev.size()) = ev;
	}
	//--------------------------------------------------
	void EmatrixToVectorDouble(Eigen::MatrixXd& emat ,Eigen::VectorXi& ev)
	{
		new (&ev) Eigen::Map<Eigen::VectorXd>(emat.data(), emat.cols()*emat.rows());
	}

	//--------------------------------------------------

	void parse2DCsvFile(std::string inputFileName, std::vector<std::vector<double>>& data) {

	    std::ifstream inputFile(inputFileName);
	    int l = 0;

	    while (inputFile) {
	        l++;
	        std::string s;
	        if (!getline(inputFile, s)) break;
	        if (s[0] != '#') {
	            std::istringstream ss(s);
	            std::vector<double> record;

	            while (ss) {
	                std::string line;
	                if (!getline(ss, line, ','))
	                    break;
	                try {
	                    record.push_back(stof(line));
	                }
	                catch (const std::invalid_argument e) {
	                    std::cout << "NaN found in file " << inputFileName << " line " << l
	                         << std::endl;
	                    e.what();
	                }
	            }

	            data.push_back(record);
	        }
	    }

	    if (!inputFile.eof()) {
	        std::cerr << "Could not read file " << inputFileName << "\n";
	        std::__throw_invalid_argument("File not found.");
	    }

	}

	void readCSV(std::string file, int rows, int cols, Eigen::MatrixXd& data)
	{
		data.resize(rows,cols);
		std::vector<std::vector<double>> data_vec;
		parse2DCsvFile(file, data_vec);

	   for(size_t i = 0; (i < rows); i++)
	   { 
	      for(size_t j = 0; (j < cols); j++)
	      {
	          data(i,j) = data_vec[i][j];
	      } 
	   }
	   data_vec.clear();
	}

	int csvWrite(const Eigen::MatrixXd& inputMatrix, const std::string& fileName, const std::streamsize dPrec)
	{
		int i, j;
		std::ofstream outputData;
		outputData.open(fileName);
		if (!outputData)
			return -1;
		outputData.precision(dPrec);
		for (i = 0; i < inputMatrix.rows(); i++) {
			for (j = 0; j < inputMatrix.cols(); j++) {
				outputData << inputMatrix(i, j);
				if (j < (inputMatrix.cols() - 1))
					outputData << ",";
			}
			if (i < (inputMatrix.rows() - 1))
				outputData << endl;
		}
		outputData.close();
		if (!outputData)
			return -1;
		return 0;
	}

	int csvWrite(double* data, int n, const std::string& fileName, const std::streamsize dPrec)
	{
		int i, j;
		std::ofstream outputData;
		outputData.open(fileName);
		if (!outputData)
			return -1;
		outputData.precision(dPrec);
		for (i = 0; i < n; i++) {
			outputData << data[i];
			if (i < (n - 1)){
				outputData << "\n";
			}else{
				outputData << endl;
			}		
		}
		outputData.close();
		if (!outputData)
			return -1;
		return 0;
	}

	void moveItemToBack(index_2d_v& v2d, size_t itemIndex)
	{
	    auto it = v2d.begin() + itemIndex;
	    std::rotate(it, it + 1, v2d.end());
	}

	void moveItemToBack(index_v& v1d, size_t itemIndex)
	{
	    auto it = v1d.begin() + itemIndex;
	    std::rotate(it, it + 1, v1d.end());
	}

	void EMatToMap(Eigen::MatrixXd& EM, bool birth, size_t bdId, index_v& ordered_nids, std::map<std::tuple<size_t, size_t>, double>& mapA)
	{
		size_t n = EM.cols();

		if (birth)
		{
			ordered_nids.pop_back();
			ordered_nids.push_back(bdId*2);
			ordered_nids.push_back(bdId*2+1);

		}else{
			ordered_nids.pop_back();
			ordered_nids.pop_back();
			ordered_nids.push_back(bdId);
		}

		for (size_t i = 0; i < n; ++i)
		{
			for (size_t j = 0; j < n; ++j)
			{
    			mapA.insert(std::make_pair(std::make_tuple(ordered_nids[i], ordered_nids[j]), EM(i,j)));
			}
		}


	}

	//--------------------------------------------------
	//Split the data ids for the bottom node birth
	void splitBnodeXids(index_v& xids, size_t v, size_t c, xinfo& xi, dinfo& di, index_2d_v& new_split_order)
	{
		// di, jth var of ith obs is *(x + p*i+j)
		new_split_order.clear();
		double *xx;//current x
		new_split_order.resize(2);
		for (auto i = xids.begin(); i != xids.end(); i++){
			xx = di.x + (*i)*di.p;
			if(xx[v] < xi[v][c]) {
				new_split_order[0].push_back(*i);
			} else {
				new_split_order[1].push_back(*i);
			}
		}
	}
	//--------------------------------------------------
	// fit function for reordering
	void fit(tree& t, xinfo& xi, size_t p, size_t n, double *x,  double* fv, index_v& xid_nid)
	{
	   xid_nid.clear();
	   tree::tree_p bn;
	   for(size_t i=0;i<n;i++) {
	      bn = t.bn(x+i*p,xi);
	      fv[i] = bn->gettheta();
		//std::cout << "in dummyutilities::fit(): " << std::endl;
		//std::cout << bn->nid() << std::endl;
	      xid_nid.push_back(bn->nid()); // vector stores corresponding bottom node nid for i-th observation.
	   }
	   	//std::cout << std::endl;
	}
	//--------------------------------------------------
	//Get ordered (ascending by nid) nid vector
	void getOrderedNid(tree::npv& bns, index_v& ordered_nids)
	{
		ordered_nids.clear();
		size_t n_nid = bns.size();
		for (size_t i = 0; i < n_nid; i++)
		{
			ordered_nids.push_back(bns[i]->nid());
		}
		std::sort (ordered_nids.begin(), ordered_nids.end()); // sort bottom nodes by ascending order.
	}
	//--------------------------------------------------
	//Reordering the data by tree fitting
	// ordered_nids: ordered tree bottom nodes by their nid
	// xid_nid: the bottom nid corresponds to the natrual order of the data
	// new_order: 2d vector, 1st dimension: ordered bottom nids. 2nd dimension: natrual order of the data.
	void reorder(index_v& xid_nid, index_v& ordered_nids, index_2d_v& new_order)
	{
	  new_order.clear();
	  size_t n_nid = ordered_nids.size();
	  new_order.resize(n_nid);
	  size_t n = xid_nid.size();
	  for (size_t i = 0; i < n; i++)
	  {
	    for (size_t j = 0; j < n_nid; j++)
	    {
	      if (xid_nid[i]==ordered_nids[j])
	      {
	        // new_order is a 2 dimensons vector, new_order[j][i] j is the order of bottom node, 
	      	//i is the orignal order of observation (data)
	        new_order[j].push_back(i); 
	        break;
	      }
	    }
	  }
	}
	//--------------------------------------------------
	//Get matrix A
	void getMatA(pinfo& pi, index_2d_v& new_order, const Eigen::MatrixXd& SigmaInv, Eigen::MatrixXd& A)
	{

		size_t n_nodes = new_order.size();
		//A.resize(n_nodes,n_nodes);
		Eigen::MatrixXd T = Eigen::MatrixXd::Zero(n_nodes,n_nodes);
		//Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n_nodes,n_nodes);

		// set prior tau
		Eigen::VectorXd t2_ev = Eigen::VectorXd::Zero(n_nodes);
		for(size_t i=0;i<n_nodes;i++) t2_ev(i)= 1/(pi.tau*pi.tau); 
		// create Q, here Q is a diagonal matrix. But generally, it can be any prior matrix.
		Eigen::MatrixXd Q = t2_ev.asDiagonal();

		Eigen::VectorXi evc;
		Eigen::VectorXi evr;
		Eigen::MatrixXd B;

		for (int i = 0; i < n_nodes; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				evc.resize(0);
				evr.resize(0);
				B.resize(0,0);
				VectorToEvector(new_order[i], evr);
				VectorToEvector(new_order[j], evc);
				igl::slice(SigmaInv,evr,evc,B);
				T(i,j)=B.sum();
			}
		}
		//
   		Eigen::MatrixXd D = (T.diagonal()).asDiagonal();
		//std::cout << "GetA T:\n" << T << std::endl;
		//std::cout << "GetA D:\n" << D << std::endl;
		A = T + T.transpose() + Q - D;

		evc.resize(0);
		evr.resize(0);
		B.resize(0,0);
		T.resize(0,0);
		D.resize(0,0);
		Q.resize(0,0);
	}
	//Get matrix A block
	void getMatA(index_2d_v& new_order, const Eigen::MatrixXd& SigmaInv, Eigen::MatrixXd& Ab, index_2d_v& block_order)
	{
		size_t n1 = new_order.size();
		size_t n2 = block_order.size();
		size_t n11;

		if (n2 > 1) // Birth
		{
			n11 = n1-1;

		}else{ //Death

			n11 = n1-2;
			
		}

		Ab.resize(n11,n2);
		Eigen::VectorXi evc;
		Eigen::VectorXi evr;
		Eigen::MatrixXd T;
		for (int i = 0; i < n11; i++)
		{
			for (int j = 0; j <n2; j++)
			{
				evc.resize(0);
				evr.resize(0);
				T.resize(0,0);
				VectorToEvector(new_order[i], evr);
				VectorToEvector(block_order[j], evc);
				igl::slice(SigmaInv,evr,evc,T);
				Ab(i,j)=T.sum();
			}
		}
		evc.resize(0);
		evr.resize(0);
		T.resize(0,0);
	}
	//--------------------------------------------------
	// Construct A from eA (the map of A entries)
	// bd_node_order nid
	void getMatAfromeA(std::map<std::tuple<size_t, size_t>, double>& eA, index_v& ordered_nids, Eigen::MatrixXd& A)
	{
		size_t n = ordered_nids.size();

		// Now you have the right matrix order. you need to construct matrix A from eA
		for (size_t i = 0; i < n; ++i)
		{
			for (size_t j = 0; j < n; ++j)
			{
				auto search = eA.find(std::make_tuple(ordered_nids[i], ordered_nids[j]));
				if (search != eA.end()) {
				    A(i,j) = search->second;
				} else {
					//printVec(ordered_nids);
				    std::cout << "Woops, not found A entry in the A map!\n";
				    std::cout << ordered_nids[i]<< "," << ordered_nids[j] << std::endl;
				}
			}
		}

	}

	//--------------------------------------------------
	//Get scalar u
	double getScalarU(index_2d_v& new_order, const Eigen::MatrixXd& SigmaInv, const Eigen::MatrixXd& B,  double *r)
	{
		//std::vector<double> vec_y(di.y,di.y+di.n);
		//Eigen::VectorXd evec_y, obs_omega;
		int n_x = SigmaInv.rows();
		Eigen::VectorXd r_ev = Eigen::VectorXd::Zero(n_x);  //Eigen vector stores the r for 
		// initialize r_ev with r
		for (int i = 0; i < n_x; ++i)
		{
			r_ev[i]=r[i];
		}
		//VectorToEvector(vec_y, evec_y); // Convert std::vector to Eigen vector
		Eigen::VectorXd obs_omega;  //Eigen vector stores the r for 
		obs_omega = SigmaInv*r_ev; // Get omega with natrual order

		size_t n_nodes = new_order.size();
		Eigen::VectorXi ev_ord; // The order of observations in each bottom node
		double u1 = 0,u2 = 0;

		Eigen::VectorXd ev_omg_i;
		Eigen::VectorXd ev_omg_j;
		for (size_t i = 0; i < n_nodes; i++)
		{
			ev_omg_i.resize(0);
			VectorToEvector(new_order[i], ev_ord); // Convert std::vector to Eigen vector
			igl::slice(obs_omega,ev_ord,ev_omg_i);

			for (size_t j = 0; j <= i; j++)
			{
				ev_omg_j.resize(0);
				VectorToEvector(new_order[j], ev_ord); // Convert std::vector to Eigen vector
				igl::slice(obs_omega,ev_ord,ev_omg_j);

				if (j<i){

					u2 = u2 + ev_omg_i.sum()*ev_omg_j.sum()*B(i,j);
					
				}else{

					u1 = u1 + ev_omg_i.sum()*ev_omg_j.sum()*B(i,j);
				}
			}
		}

		ev_omg_i.resize(0);
		ev_omg_j.resize(0);
		obs_omega.resize(0);
		ev_ord.resize(0);
		//std::cout << "\nScalar U:" << 2*u2 + u1 << "\n" << std::endl;

		return 2*u2 + u1;
	}
	//--------------------------------------------------
	//Get birth marginal likelihood ratio given T
	double getLogMglrBirth(tree& x, pinfo& pi, double *r, const Eigen::MatrixXd& SigmaInv, index_2d_v& new_order, index_2d_v& new_split_order, size_t birth_node_order, index_v& ordered_nids, bool sigmaNoChange, Eigen::MatrixXd& mA_i1)
	{	
		size_t n_nodes = new_order.size();
		Eigen::MatrixXd A_i(n_nodes,n_nodes), invA_i(n_nodes,n_nodes);
		// Get the inverse of A^i: invA_i
		// put the birth node (given by birth_node_order) to the last one
		moveItemToBack(new_order, birth_node_order);
		moveItemToBack(ordered_nids, birth_node_order);

		// if Sigma is constant
		if (sigmaNoChange)
		{
			std::map<std::tuple<size_t, size_t>, double> eA_i;
			x.geteA(eA_i);
			getMatAfromeA(eA_i, ordered_nids, A_i);
			invA_i = A_i.inverse();

		}else{ // if Sigam is changing in every iteration, we have to calculate A and A inverse every time

			getMatA(pi, new_order, SigmaInv, A_i);
			invA_i = A_i.inverse();

		}

		// Create inverse of A^i extension: invA_i_ex
		Eigen::MatrixXd invA_i_ex(n_nodes+1,n_nodes+1);
		invA_i_ex << invA_i, invA_i.rightCols(1), invA_i.bottomRows(1), invA_i(n_nodes-1,n_nodes-1);
		//std::cout << "invA_i_ex:\n" << invA_i_ex << std::endl;
	
		// Get the A^(i+1)for inverse
		//Get A_i1 by block matrices
		Eigen::MatrixXd A_i1(n_nodes+1,n_nodes+1);
		Eigen::MatrixXd Ab1,Ab2;
		getMatA(new_order, SigmaInv, Ab1, new_split_order);
		getMatA(pi, new_split_order, SigmaInv, Ab2);

		A_i1 << A_i.block(0,0,n_nodes-1,n_nodes-1), Ab1, Ab1.transpose(), Ab2;  // we can not updata x.A here, because we don't know the tree will birth or not.
		mA_i1 = A_i1;

		// Create matrix B = invA_i1 - invA_i_ex
		Eigen::MatrixXd B = A_i1.inverse() - invA_i_ex;
		//std::cout << "B:\n" << B << "\n" << std::endl;

		// Calculate u
		double u = getScalarU(new_order, SigmaInv, B, r);
		// Calculate the marginal likelihood ritao
		//std::cout << "A_i determinant:" << A_i.determinant() << "       A_i1 determinant:" << A_i1.determinant() << "\n" << std::endl;
		double aid, ai1d;
		aid = A_i.determinant();
		ai1d = A_i1.determinant();

		//clean
		invA_i.resize(0,0);
		Ab1.resize(0,0);
		Ab2.resize(0,0);
		invA_i_ex.resize(0,0);
		B.resize(0,0);
		A_i.resize(0,0);
		A_i1.resize(0,0);
		
		// Calculate the marginal likelihood ritao
		if (aid>0 && ai1d>0)
		{	//Because Q is a diagonal matrix so we have -log(pi.tau)
			return -log(pi.tau) + 0.5*(log(aid)- log(ai1d)) + u/2;

		}else{

			return -std::numeric_limits<double>::infinity();
		}
	}
	//--------------------------------------------------
	//Get death marginal likelihood ratio given T
	double getLogMglrDeath(tree& x, pinfo& pi, double *r, const Eigen::MatrixXd& SigmaInv, index_2d_v& new_order, size_t death_nodes_order, index_v& ordered_nids, bool sigmaNoChange, Eigen::MatrixXd& mA_i1)
	{
		size_t n_nodes = new_order.size();
		// Get the A^i for inverse
		for (size_t i = 0; i < 2; i++) //we need to move two nodes to the back
		{
			// this is because after first rotate, the order of the second item reduces 1
			moveItemToBack(new_order, death_nodes_order);
			moveItemToBack(ordered_nids, death_nodes_order);
		}
		//std::cout << "I'm here!!!! ========= 2" << std::endl;
		//printVec(ordered_nids);
		//std::cout << "I'm here!!!! ========= 2" << std::endl;
		Eigen::MatrixXd A_i(n_nodes,n_nodes), invA_i(n_nodes,n_nodes);
		// if Sigma is constant
		if (sigmaNoChange)
		{
			std::map<std::tuple<size_t, size_t>, double> eA_i;
			x.geteA(eA_i);
			getMatAfromeA(eA_i, ordered_nids, A_i);
			invA_i = A_i.inverse();

		}else{ // if Sigam is changing in every iteration, we have to calculate A and A inverse every time
			getMatA(pi, new_order, SigmaInv, A_i);
			invA_i = A_i.inverse();
		}

		// Get the inverse of A^(i+1): invA_i1
		Eigen::MatrixXd A_i1(n_nodes-1,n_nodes-1),invA_i1(n_nodes-1,n_nodes-1);

		index_v end, pre_end;
		pre_end = new_order[n_nodes-2];
		end = new_order[n_nodes-1];
		std::move(end.begin(),end.end(),std::back_inserter(pre_end));
		//std::cout << "\n# of new_order:" << n_nodes << "    after move:    " << new_order.size() << std::endl;

		index_2d_v concatenate_order;
		concatenate_order.push_back(pre_end);
		Eigen::MatrixXd Ab1,Ab2;
		getMatA(new_order, SigmaInv, Ab1, concatenate_order);
		getMatA(pi, concatenate_order, SigmaInv, Ab2);
		A_i1 << A_i.block(0,0,n_nodes-2,n_nodes-2), Ab1, Ab1.transpose(), Ab2;
		mA_i1 = A_i1;

		//clean
		Ab1.resize(0,0);
		Ab2.resize(0,0);

		invA_i1 = A_i1.inverse();

		// Create inverse of A^(i+1) extension: invA_i1_ex
		Eigen::MatrixXd invA_i1_ex(n_nodes,n_nodes);
		invA_i1_ex << invA_i1, invA_i1.rightCols(1), invA_i1.bottomRows(1), invA_i1(n_nodes-2,n_nodes-2);

		invA_i1.resize(0,0);

		// Create matrix B = invA_i1_ex - invA_i
		Eigen::MatrixXd B = invA_i1_ex - invA_i;

		invA_i1_ex.resize(0,0);

		// Calculate u
		double u = getScalarU(new_order, SigmaInv, B, r);

		B.resize(0,0);

		double aid, ai1d;
		aid = A_i.determinant();
		ai1d = A_i1.determinant();

		A_i.resize(0,0);
		A_i1.resize(0,0);

		// Calculate the marginal likelihood ritao
		if (aid>0 && ai1d>0)
		{	//Because Q is a diagonal matrix so we have log(pi.tau)
			return log(pi.tau) + 0.5*(log(aid)- log(ai1d)) + u/2;

		}else{

			return -std::numeric_limits<double>::infinity();
		}

	}
	//--------------------------------------------------
	//Metropolis-Hasting step for birth and death
	bool bd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double *r, const Eigen::MatrixXd& SigmaInv,
		std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen, index_v& xid_nid, bool sigmaNoChange)
	{
		//std::cout << "enter new bd!===================" << std::endl;
		tree::npv goodbots;  //nodes we could birth at (split on)
		double PBx = getpb(x,xi,pi,goodbots); //prob of a birth at x
		//--------------------------------------------------
		//get new order
		index_v ordered_nids;
		index_2d_v new_order;
		tree::npv bnv; //all the bottom nodes
		x.getbots(bnv);
		getOrderedNid(bnv, ordered_nids);
		bnv.clear();

		Eigen::MatrixXd matA_i1;

		//if((gen.uniform() < PBx)&&(treesize < 5)) { //do birth or death
		if(gen.uniform() < PBx) { //do birth or death

			reorder(xid_nid, ordered_nids, new_order);
			//--------------------------------------------------
			//draw proposal
			tree::tree_p nx; //bottom node
			size_t v,c; //variable and cutpoint
			double pr; //part of metropolis ratio from proposal and prior
			bprop(x,xi,pi,goodbots,PBx,nx,v,c,pr,nv,pv,aug,gen);
			size_t bnid_birth = nx->nid();
			//std::cout << "bnid_birth:" << bnid_birth << std::endl;
			//printVec(ordered_nids);
			//--------------------------------------------------
			index_v xids;
			index_2d_v new_split_order;
			size_t birth_node_order; // the order of birth node in all the ordered bottom nodes
			size_t n_onid = ordered_nids.size();
			for (size_t i = 0; i < n_onid; i++)
			{
				if (bnid_birth==ordered_nids[i])
				{
					xids = new_order[i];
					birth_node_order = i;

				}
			}
			splitBnodeXids(xids, v, c, xi, di, new_split_order);
			size_t nr,nl; //counts in proposed bots
			nl = new_split_order[0].size();
			nr = new_split_order[1].size();
			//compute alpha
			double alpha=0.0, lalpha=0.0;

			if((nl>=5) && (nr>=5)) { //cludge?
				double logmglrb = getLogMglrBirth(x, pi, r, SigmaInv, new_order, new_split_order, birth_node_order, ordered_nids, sigmaNoChange, matA_i1);
				alpha=1.0;
				lalpha = log(pr) + logmglrb;
				//std::cout << "lalpha (birth):" << lalpha << std::endl;
				lalpha = std::min(0.0,lalpha);
			}

			//--------------------------------------------------
			//try metrop
			double uu = gen.uniform();
			bool dostep = (alpha > 0) && (log(uu) < lalpha);
			if(dostep) {
				// set the ltheta and rtheta both equal to 0, it's because no matter what value you set
				// you always draw theta for all bottom nodes in the following step drmu()
				x.birthp(nx,v,c,0,0);
				nv[v]++;
				// here, we need to update eA for birth
				std::map<std::tuple<size_t, size_t>, double> eA_i1;
				EMatToMap(matA_i1, true, bnid_birth, ordered_nids ,eA_i1);
				x.reseteA(eA_i1);
		    	return true;
			} else {
		    	return false;
			}
		} else {

			// do death
			reorder(xid_nid, ordered_nids, new_order);
			//--------------------------------------------------
			//draw proposal
			double pr;  //part of metropolis ratio from proposal and prior
			tree::tree_p nx; //nog node to death at
			dprop(x,xi,pi,goodbots,PBx,nx,pr,gen);

			size_t n_onid = ordered_nids.size();
			size_t deathId = nx->nid();
			size_t deathFirstChildNid = 2*deathId; // the order of first death node in all the ordered bottom nodes, the next one must be death
			size_t death_nodes_order;

			for (size_t i = 0; i < n_onid; i++)
			{
				if (deathFirstChildNid==ordered_nids[i])
				{
					death_nodes_order = i;
				}
			}

			//--------------------------------------------------
			// calculate lalpha
			double logmglrd = getLogMglrDeath(x, pi, r, SigmaInv, new_order, death_nodes_order, ordered_nids, sigmaNoChange, matA_i1);
			double lalpha = log(pr) + logmglrd;
			//std::cout << "lalpha (death):" << lalpha << std::endl;
			lalpha = std::min(0.0,lalpha);
			//std::cout << "lalpha (death):" << lalpha << std::endl;

			//--------------------------------------------------
			//try metrop
			double death_gen = log(gen.uniform());
			//std::cout << "death_gen:" << death_gen << std::endl;
			if( death_gen < lalpha) {
				nv[nx->getv()]--;
				x.deathp(nx,0);
				// here, we need to update eA for death
				std::map<std::tuple<size_t, size_t>, double> eA_i1;
				EMatToMap(matA_i1, false, deathId, ordered_nids, eA_i1);
				x.reseteA(eA_i1);
				return true;
			} else {
				return false;
			}
		}
	}
	//--------------------------------------------------
	//draw all the dependent bottom node mu's
	void drmu(tree& t, pinfo& pi, double *r, const Eigen::MatrixXd& SigmaInv, index_v& xid_nid ,bool sigmaNoChange)
	{
		//std::cout << "\n\nNow draw mus:\n" << std::endl;

		//std::cout << "print xid_nid:\n" << std::endl;
		//printVec(xid_nid);

		//t.pr();

		tree::npv bnv; //all the bottom nodes
		t.getbots(bnv);
		//--------------------------------------------------
		//get new order
		index_v ordered_nids;
		index_2d_v new_order;
		getOrderedNid(bnv, ordered_nids);

		size_t n_btn = ordered_nids.size();
		//std::cout << "print bottom nodes by ordered nid:\n" << std::endl;
		//printVec(ordered_nids);

		reorder(xid_nid, ordered_nids, new_order);
		//std::cout << "print new data order in bottom nodes:\n" << std::endl;
		//printVec(ordered_nids);
		//--------------------------------------------------
		// Calculation
		//std::cout << "\nSigmaInv:\n" << SigmaInv << std::endl;
		//std::cout << "\nNow calculate posterior inverse covariance matrix:\n" << std::endl;
		Eigen::MatrixXd A(n_btn,n_btn), InvA;

		if (sigmaNoChange)
		{
			std::map<std::tuple<size_t, size_t>, double> eA;
			t.geteA(eA);
			getMatAfromeA(eA, ordered_nids, A);
			//getMatA(pi, new_order, SigmaInv, A2);// get the inverse covariance matrix
			//std::cout << "A map method:===================\n" << A << std::endl;
			//std::cout << "Old method:===================\n" << A2 << std::endl;

		}else{
			getMatA(pi, new_order, SigmaInv, A);// get the inverse covariance matrix
			//std::cout << "I'm here!!!! ========= 4\n" << SigmaInv.sum() << std::endl;
		}

	    //std::map<std::tuple<size_t, size_t>, double>::iterator it1 = eA.begin();
	    //while(it1 != eA.end())
	    //{
	        //std::cout<< std::get<0>(it1->first) << ", " << std::get<1>(it1->first) <<" :: "<< it1->second<<std::endl;
	        //it1++;
	    //}

		//std::cout << "I'm here!!!! ========= 5\n" << A << std::endl;

		InvA = A.inverse();
		//std::cout << "I'm here!!!! ========= 6\n" << InvA << std::endl;
		//exit(0);
		//clean
		A.resize(0,0);

		// calculate D_P * Sigma_P * R_P
		size_t n_x = xid_nid.size();
		Eigen::VectorXd r_ev = Eigen::VectorXd::Zero(n_x);  //Eigen vector stores the r for 
		for (int i = 0; i < n_x; ++i)
		{
			r_ev[i]=r[i];
		}

		Eigen::VectorXd obs_omega;

		obs_omega = SigmaInv*r_ev; // Get omega with natrual order
		//std::cout << "obs_omega:\n" << obs_omega << std::endl;

		size_t n_nodes = new_order.size();
		//Eigen::VectorXd temp(n_nodes), avgtemp(n_nodes),avgy(n_nodes);
		Eigen::VectorXd temp(n_nodes);
		Eigen::VectorXi num_nodes(n_nodes);

		Eigen::VectorXi ev_ord;
		Eigen::VectorXd ev_omg_i,evec_y_i;

		for (size_t i = 0; i < n_nodes; i++)
		{
			ev_ord.resize(0);
			ev_omg_i.resize(0);
			VectorToEvector(new_order[i], ev_ord); // Convert std::vector to Eigen vector
			igl::slice(obs_omega,ev_ord,ev_omg_i);
			igl::slice(r_ev,ev_ord,evec_y_i);
			temp(i)=ev_omg_i.sum();
			//std::cout << "ev_omg_i:\n" << ev_omg_i << std::endl;
		}
		//clean
		ev_ord.resize(0);
		ev_omg_i.resize(0);

		Eigen::VectorXd mu;
		mu = InvA*temp;
		//std::cout << "mu===============:\n" << mu << std::endl;

		// inintial posterior
		Mvn genMvn(mu,InvA); 
		Eigen::MatrixXd mus;

		//draw sample
		genMvn.sample(mus, 1);
		//std::cout << "drawmu===================:\n" << mus << std::endl;
		for (size_t i = 0; i < n_nodes; i++)
		{
			(t.getptr(ordered_nids[i]))->settheta(mus(i));
		}
	}

}//namespace dummyutilities