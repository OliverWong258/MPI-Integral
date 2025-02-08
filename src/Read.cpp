#include"Read.h"

void MyRead::ReadInput(string pathname){
    arg2meth["isHexahedral"] = "0";
    arg2meth["lx"] = "1000";
    arg2meth["ly"] = "1000";
    arg2meth["lz"] = "1000";
    arg2meth["thetaxy"] = "0";
    arg2meth["thetayz"] = "0";
    arg2meth["thetazx"] = "0";
    arg2meth["support_SH"] = "0";
    arg2meth["diago_lib"] = "lapack";
    arg2meth["support_Periodic_Boundary"] = "0";
    arg2meth["multi_parallel_strategies"] = "0";
    
    ifstream ifs(pathname);
    if (!ifs) {//fail to open the input file
		cout << "can not open the Input file" << endl;
		cout << pathname << endl;
		exit(1);
	}
    string str;
	string argname = "", methname = "";
    while (getline(ifs, str)) {//read in the argument names and method names
		istringstream is(str);
		is >> argname >> methname;
		if(argname[0] == '#')//only commentary
			continue;
		arg2meth[argname] = methname;
	}
	lx = stod(arg2meth["lx"]);
	ly = stod(arg2meth["ly"]);
	lz = stod(arg2meth["lz"]);
	ifs.close();
}

void MyRead::ReadPoints(){
    ifstream ifs(arg2meth["points_path"]);
    if (!ifs) {//fail to open the file
		cout << "can not open the Points file" << endl;
		cout << arg2meth["points_path"] << endl;
		exit(1);
	}
    string str;
	string x1, x2, x3;
	point_vec = new point[50];
    while (getline(ifs, str)) {//read in the argument names and method names
        replace(str.begin(), str.end(), ',', ' ');
		replace(str.begin(), str.end(), '(', ' ');
		replace(str.begin(), str.end(), ')', ' ');
		istringstream is(str);
		is >> x1 >> x2 >> x3;
		point tmp;
		tmp.x = stod(x1);
		tmp.y = stod(x2);
		tmp.z = stod(x3);
		point_vec[point_num] = tmp;
		++point_num;
	}
	ifs.close();
}

void MyRead::ReadV(){
	ifstream ifs(arg2meth["v_path"]);
	if (!ifs) {//fail to open the file
		cout << "can not open the V file" << endl;
		cout << arg2meth["V_path"] << endl;
		exit(1);
	}
    string str;
	string str1, str2;
	for(int i = 0;i < 3;++i){
		getline(ifs, str);
		istringstream is(str);
		is >> str1 >> str2;
		if(str1 == "nx"){
			nx = stoi(str2);
		}
		else if(str1 == "ny"){
			ny = stoi(str2);
		}
		else{
			nz = stoi(str2);
		}
	}
	getline(ifs, str);
	long long num = nx*ny*nz;
	double tmp;
	V = new double[num];
	for(long long i = 0;i < num;++i){
		ifs >> tmp;
		V[i] = tmp;
	}
	ifs.close();
}

void MyRead::ReadDistr(){
	ifstream ifs(arg2meth["distribution_path"]);
	if (!ifs) {//fail to open the file
		cout << "can not open the distribution file" << endl;
		cout << arg2meth["distribution_path"] << endl;
		exit(1);
	}
	string str1, str2, str;
	for(int i = 0;i < 4;++i){
		ifs >> str1 >> str2;
		if(str1 == "cutoff")
			cutoff = stod(str2);
		else if(str1 == "dr")
			dr = stod(str2);
		else if(str1 == "mesh")
			mesh = stoi(str2);
		else if(str1 == "l")
			l = stoi(str2);
	}
	ifs >> str;
	str = "";
	double tmp;
	distr = new double[mesh];
	int index = 0;
	while(getline(ifs, str)){
		replace(str.begin(), str.end(), ',', ' ');
		istringstream is(str);
		while(is >> tmp){
			distr[index] = tmp;
			index++;
		}
	}
	ifs.close();
}