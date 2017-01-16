class masspoint
{
   public:
      double getMS( void );
      double getMZ( void );
      std::vector<double > getPos( void );
      std::vector<double > getChi( void );
      masspoint(double  ms, double mz);  // This is the constructor
      int marginalize(std::vector<double > *bfChi, std::vector<double > *chi);
   private:
      double Ms;
      double Mz;
      std::vector<std::vector<double >  > minvalues;
      int numsteps = 30;
      int numruns = 2;
// 3: 0.01->0.75 , 0.15->0.6
      std::vector<double > lbound = {0.65, 0.5,  0.125, 0.15, 0.4,  0.35,  0.01,  0.4, 0.5 ,0,0,0,0};
      std::vector<double > ubound = {0.9,  0.65, 0.175, 0.6 ,  0.75,  0.8,  0.7,   0.8,  0.9 ,2*3,14159,2*3.14159,2*3.14159,2*3.14159};
      std::vector<double > stepsize = {(ubound[0]-lbound[0])/numsteps,(ubound[1]-lbound[1])/numsteps,(ubound[2]-lbound[2])/numsteps,(ubound[3]-lbound[3])/numsteps,(ubound[4]-lbound[4])/numsteps,(ubound[5]-lbound[5])/numsteps,(ubound[6]-lbound[6])/numsteps,(ubound[7]-lbound[7])/numsteps,(ubound[8]-lbound[8])/numsteps,(ubound[9]-lbound[9])/numsteps,(ubound[10]-lbound[10])/numsteps,(ubound[11]-lbound[11])/numsteps,(ubound[12]-lbound[12])/numsteps};

};

masspoint::masspoint(double ms, double mz){
    Ms=ms;
    Mz=mz;
}

double masspoint::getMS( void ){    return Ms;}
double masspoint::getMZ( void ){    return Mz;}

int masspoint::marginalize(std::vector<double > *bfChi, std::vector<double > *chi){

	std::ofstream teststream("melement_taubound_r_" + std::to_string(whichele)+".dat");

	double PI =3.14159265359;
	nlopt::opt ue1(nlopt::GN_ISRES,12);
	nlopt::opt ue1Loc(nlopt::LN_COBYLA,12);

	ue1.set_maxeval(50000);
	//ue1.set_xtol_abs(1e-4);
	ue1Loc.set_xtol_abs(1e-4);

	std::vector<double > lb2 = {0,0,0,0,0,0,0,0,0,0,0,0};
	std::vector<double > ub2={1,1,1,1,1,1,1,1,2*PI,2*PI,2*PI,2*PI};
	ue1.set_lower_bounds(lb2);
	ue1.set_upper_bounds(ub2);
	ue1Loc.set_lower_bounds(lb2);
	ue1Loc.set_upper_bounds(ub2);

	my_constraint_data data[2] = { {double (whichele),0.0}, {-1,1} };
	my_func_data ddata[2]={ {double (whichele),0.0}, {-1,1} };

	ue1.set_min_objective(testfue1, &ddata[0]);
	ue1Loc.set_min_objective(testfue1, &ddata[0]);

	ue1.add_inequality_constraint(myvconstraint1WHI, &data[0], 1e-8);
	ue1.add_inequality_constraint(myvconstraint2WHI, &data[0], 1e-8);
	ue1.add_inequality_constraint(myvconstraint3WHI, &data[0], 1e-8);
	ue1.add_inequality_constraint(geoCon1variable,	 &data[0],1e-8);
	ue1.add_inequality_constraint(geoCon2variable,	 &data[0],1e-8);
	ue1.add_inequality_constraint(geoCon3variable,	 &data[0],1e-8);
	ue1.add_inequality_constraint(geoCon4variable,   &data[0],1e-8);
	ue1.add_inequality_constraint(geoCon5variable,   &data[0],1e-8);
	ue1.add_inequality_constraint(geoCon6variable,   &data[0],1e-8);
	ue1.add_inequality_constraint(minEnormWHI,       &data[0],1e-8);

	ue1Loc.add_inequality_constraint(myvconstraint1WHI,	 &data[0], 1e-8);
	ue1Loc.add_inequality_constraint(myvconstraint2WHI,	 &data[0], 1e-8);
	ue1Loc.add_inequality_constraint(myvconstraint3WHI,	 &data[0], 1e-8);
	ue1Loc.add_inequality_constraint(geoCon1variable,	&data[0],  1e-8);
	ue1Loc.add_inequality_constraint(geoCon2variable,	&data[0],  1e-8);
	ue1Loc.add_inequality_constraint(geoCon3variable,	&data[0],  1e-8);
	ue1Loc.add_inequality_constraint(geoCon4variable,	&data[0],  1e-8);
	ue1Loc.add_inequality_constraint(geoCon5variable,	&data[0],  1e-8);
	ue1Loc.add_inequality_constraint(geoCon6variable,	&data[0],  1e-8);
	ue1Loc.add_inequality_constraint(minEnormWHI,		&data[0],  1e-8);

	std::vector<double> xe1= {0.8235,0.546549,0.151778,0.496027, 0.56412,0.660095,0.275153,0.618913,0.735688,PI,1.99*PI,PI,1.99*PI};
	//std::vector<double > xe1 = {0.823558, 0.546549, 0.151778,0.450511, 0.467655, 0.760486, 0.344663, 0.694682, 0.631367,PI,1.99*PI,PI,1.99*PI};
	xe1.erase (xe1.begin()+whichele);
	std::vector<double > xalways = xe1;
	
	std::vector<nlopt::result > results;
	results.reserve(numruns);

	std::vector<std::vector<double > > xtmp;
	for(int j=0;j<numruns;j++){xtmp.push_back(xe1); };
	std::vector<std::vector<double > > xtmpLoc;
	for(int j=0;j<numruns;j++){xtmpLoc.push_back(xe1); };
	std::vector<std::vector<double > > xtmpLocSole;
	for(int j=0;j<numruns;j++){xtmpLocSole.push_back(xe1); };


	std::vector<double >  min(numruns);
	std::vector<double >  minLoc(numruns);
	std::vector<double >  minLocSole(numruns);


	for(int i=0; ubound[whichele]-stepsize[whichele]*i >= lbound[whichele] ;i++){
		double k =  ubound[whichele]-stepsize[whichele]*i ;//lbound[whichele]+stepsize[whichele]*i;// ubound[whichele]-stepsize[whichele]*i ;
	 	ddata[0] =  {double (whichele),k};
	 	data[0] =   {double (whichele),k};
		

	//std::vector<double> xe1= {0.8235,0.546549,0.151778,0.496027, 0.56412,0.660095,0.275153,0.618913,0.735688,PI,1.99*PI,PI,1.99*PI};
	//xe1.erase (xe1.begin()+whichele); // Delete if you want to store anything good.


		for(int j=0;j<numruns;j++){ //reinitialise all "numrun" min vectors with best from last time.
			xtmp[j]=xe1;
			xtmpLoc[j]=xe1;
			xtmpLocSole[j]=xe1;
		};


		for(int j=0;j<numruns;j++){
			nlopt::srand_time();
			results[j]  = ue1.optimize(xtmp[j], min[j]);	
			nlopt::srand_time();	
			results[j]  = ue1Loc.optimize(xtmpLoc[j], minLoc[j]);	
			nlopt::srand_time();
			results[j]  = ue1Loc.optimize(xtmpLocSole[j], minLocSole[j]);		
		//	xtmpLocSole[j]=xalways;
		//	results[j]  = ue1Loc.optimize(xtmpLocSole[j], minLocSole[j]);	
		};


		double minnextGlob = *std::min_element(min.begin(), min.end());
		double minnextLoc  = *std::min_element(minLoc.begin(), minLoc.end());
		double minnextLocSole  = *std::min_element(minLocSole.begin(), minLocSole.end());
		double minnext_temp =std::min(minnextGlob,minnextLoc);
		double minnext =   std::min(minnext_temp,minnextLocSole);

		for(int j=0;j<numruns;j++){
			if (minnext == min[j]){
				xe1 = xtmp[j];
				std::cout<<"Glob"<<std::endl;
			} else if (minnext == minLoc[j]){
				xe1 = xtmpLoc[j];
				std::cout<<"Local"<<std::endl;
			} else if (minnext == minLocSole[j]){
				xe1 = xtmpLocSole[j];
				std::cout<<"Sole"<<std::endl;
			}
		};

		position.push_back(k);
		chi2.push_back(minnext);
		
		minvalues.push_back( insertval(xe1, whichele, k)); // insert best fit into the dataset

	    	//result2 = ue1Loc.optimize(xtmp, minue2Loc);
		//result3 = ue1Loc.optimize(xe1, minue1Loc);

		//result3b = ue1.optimize(xe1, minue1);
		//result5b = ue1Loc.optimize(xe1, minue1Loc);
	    	
	   // std::cout<<minue1Loc<<" "<<xe1[0]<<" "<<xe1[1]<<" "<<k<<" "<<xe1[2]<<" "<<xe1[3]<<" "<<xe1[4]<<" "<<xe1[5]<<" "<<xe1[6]<<" "<<xe1[7]<<" "<<xe1[8]<<" "<<xe1[9]<<" "<<xe1[10]<<" "<<xe1[11]<<std::endl;
		//dispvec(min);
		std::vector<double > ye1 = insertval(xe1, whichele, k);
		std::cout<<"which: "<<whichele<<" Pos: "<<k<<" Min: "<<minnext<<std::endl;
		teststream<<k<<" "<<minnext<<" ";
		teststream<<ye1[0]<<" "<<ye1[1]<<" "<<ye1[2]<<" "<<ye1[3]<<" "<<ye1[4]<<"  "<<ye1[5]<<" "<<ye1[6]<<" "<<ye1[7]<<" "<<ye1[8]<<" "<<ye1[9]<<" "<<ye1[10]<<" "<<ye1[11]<<"  "<<ye1[12]<<std::endl;
		
	 
	}
	
	(*pos)=position;
	(*chi)=chi2;

	    return 1;
}

