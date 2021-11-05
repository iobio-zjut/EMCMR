#include "MC.h"
#include "GetVdwRadius.h"
#include <time.h>
#include "Rotbuilder.h"
#include "CaDensityx.h"
#include "randomx.h"
#include <unistd.h>

using namespace std;

/*
random_device rdTRNG;
mt19937 mtRNG;
std::uniform_real_distribution<double> fxdist;
std::uniform_int_distribution<int> ixdist; 

double randf0and1(){
    //The line below isn't necessary since the fdist distribution should always be between 0 and 1.
    //But it's there (for now) due to clarity.
    fxdist.param(std::uniform_real_distribution<double>(0.0, 1.0).param());
    return fxdist(mtRNG);
}
int randIntCustom(int a, int b){ 
    ixdist.param(std::uniform_int_distribution<int>(a, b).param());
    return ixdist(mtRNG);
}  */

std::string
name2eltxx( std::string line ) {
    std::string atmid = line.substr(12,4);
    while ( !atmid.empty() && atmid[0] == ' ' ) atmid = atmid.substr(1,atmid.size()-1);
    while ( !atmid.empty() && atmid[atmid.size()-1] == ' ' ) atmid = atmid.substr(0,atmid.size()-1);
    std::string resname = line.substr(17,3);
    while ( !resname.empty() && resname[0] == ' ' ) resname = resname.substr(1,resname.size()-1);
    while ( !resname.empty() && resname[resname.size()-1] == ' ' ) resname = resname.substr(0,resname.size()-1);

    std::string type;

    if ( line.substr(0,4) == "ATOM" ) {
        type = atmid.substr(0,2);
        if ( isdigit(type[0]) ) {
            // sometimes non-standard files have, e.g 11HH
            if ( !isdigit(type[1]) ) type = atmid.substr(1,1);
            else type = atmid.substr(2,1);
        } else if ( (line[12] == ' ' && type!="Zn" && type!="Fe" && type!="ZN" && type!="FE")
                || isdigit(type[1]) ) {
            type = atmid.substr(0,1);     // one-character element
        }

        if ( resname.substr(0,2) == "AS" || resname[0] == 'N' ) {
            if ( atmid == "AD1" ) type = "O";
            if ( atmid == "AD2" ) type = "N";
        }
        if ( resname.substr(0,3) == "HIS" || resname[0] == 'H' ) {
            if ( atmid == "AD1" || atmid == "AE2" ) type = "N";
            if ( atmid == "AE1" || atmid == "AD2" ) type = "C";
        }
        if ( resname.substr(0,2) == "GL" || resname[0] == 'Q' ) {
            if ( atmid == "AE1" ) type = "O";
            if ( atmid == "AE2" ) type = "N";
        }
        if ( atmid.substr(0,2) == "HH" ) { // ARG
            type = "H";
        }
        if ( atmid.substr(0,2) == "HD" || atmid.substr(0,2) == "HE" || atmid.substr(0,2) == "HG" ) {
            type = "H";
        }
    } else {
        if ( isalpha(atmid[0]) ) {
            if ( atmid.size() > 2 && (atmid[2] == '\0' || atmid[2] == ' ') ) {
                type = atmid.substr(0,2);
            } else if ( atmid[0] == 'A' ) { // alpha prefix
                type = atmid.substr(1, atmid.size() - 1);
            } else {
                type = atmid.substr(0,1);
            }
        } else if ( atmid[0] == ' ' ) {
            type = atmid.substr(1,1); // one char element
        } else {
            type = atmid.substr(1,2);
        }

        if ( atmid == resname ) {
            type = atmid;
            if ( type.size() == 2 ) type[1] = toupper(type[1]);
        } else if ( resname == "ADR" || resname == "COA" || resname == "FAD" ||
                resname == "GPG" || resname == "NAD" || resname == "NAL" ||
                resname == "NDP" || resname == "ABA" )  {
            if ( type.size() > 1 ) type = type.substr(0,1);
        } else if ( isdigit(type[0]) ) {
            type = type.substr(1,1);
        } else if ( type.size() > 1 && isdigit(type[1]) ) {
            type = type.substr(0,1);
        } else if ( type.size() > 1 && isalpha(type[1]) ) {
            if ( type[0] == 'O' && type[1] == 'H' ) {
                type = type.substr(0,1); // no "Oh" element (e.g. 1MBN)
            } else if ( islower(type[1]) ) {
                type[1] = toupper(type[1]);
            }
        }
    }
    return type;
}

vector<string> string_splitxx(string const & in, char splitchar)
{
    vector<string> parts;
    int i=0,j=0;
    while( j!= std::string::npos)
    {
        j=in.find(splitchar,i);
        std::string const part = in.substr(i,j-i);
        parts.push_back( part);
        i=j+1;
    }
    return parts;
}

void writePDBcoordss(std::string filename, poseCoords &atmlist)
{
    std::ifstream inpdb(filename.c_str());
//  cout<< filename.c_str()<<endl;
    std::string buf;
    vector<string> tmp_str;
    tmp_str = string_splitxx(filename,'.');
    string filenamex = tmp_str[0] + "_out.pdb";
    std::ofstream outpdb;
    outpdb.open(filenamex.c_str());
    int tmpx=0;

    while ( std::getline(inpdb, buf ) ) {
//        if ( buf.substr(0,4)!="ATOM" && buf.substr(0,6)!="HETATM" ) continue;
    	if(buf.length()<54) continue;
        if ( buf.substr(0,4) =="ATOM" && buf.substr(13,1)!="H") {

	        poseCoord atom_i;
	        string tmp_A,tmp_B;

	        tmp_A = buf.substr(0,30);
	        atom_i.x_.push_back(atof(buf.substr(30,8).c_str()));
	        atom_i.x_.push_back(atof(buf.substr(38,8).c_str()));
	        atom_i.x_.push_back(atof(buf.substr(46,8).c_str()));
	//        tmp_B = buf.substr(54,6);
	//        atom_i.B_ = atof(buf.substr(60,6).c_str());
	    //    cout<<"YYY"<<endl;


	//        atom_i.elt_ = name2elt( buf ); // horrible hacky logic mapping name->elt (could use PDB fields 76-77 if on by default)
	//        if ( atom_i.elt_ == "H" ) continue;

	        outpdb<<setw(30)<<tmp_A;
	        outpdb<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<atmlist[tmpx].x_[0];
	        outpdb<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<atmlist[tmpx].x_[1];
	        outpdb<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<atmlist[tmpx].x_[2];
	//        outpdb<<setw(6)<<tmp_B;
	//        outpdb<<setw(6)<<setiosflags(ios::fixed)<<setprecision(2)<<atmlist[tmpx].B_;
	        outpdb<<endl;
	        tmpx = tmpx + 1;        
	//        atmlist.push_back( atom_i );
    	}
    }
    outpdb.close();
}

// quick and dirty PDB read where we only care about heavyatom locations and atom ids
void
readPDBcoordss(std::string filename, poseCoords &atmlist) {
    std::ifstream inpdb(filename.c_str());
//  cout<< filename.c_str()<<endl;
    std::string buf;

    int nn=0;
    while ( !inpdb.eof() && std::getline(inpdb, buf) ) {
//    	cout<<nn << " ";
    	if(buf.length()<54) continue;
//    	if ( buf.substr(0,4) =="ATOM" && buf.substr(12,4) ==" CA " ) {
    	if ( buf.substr(0,4) =="ATOM" && buf.substr(13,1) !="H" && buf.substr(13,3) !="OXT") {
	        poseCoord atom_i;

	        atom_i.x_.push_back(atof(buf.substr(30,8).c_str()));
	        atom_i.x_.push_back(atof(buf.substr(38,8).c_str()));
	        atom_i.x_.push_back(atof(buf.substr(46,8).c_str()));
	        atom_i.B_ = 0.0;// atof(buf.substr(60,6).c_str())
	        atom_i.elt_ = buf.substr(12,4);

//	        nn=nn+1;
//	        cout<< nn<<" ";
	        atmlist.push_back( atom_i );
    	}
    }
}

vector<float> Generate_random_point(float A,float B,float r)
{
	vector<float> coordx;
	float u=A+(B-A)*rand()/(RAND_MAX+1.0) ; // [a,b] random number,(rand()%(b-a+1))+a
	float v=A+(B-A)*rand()/(RAND_MAX+1.0);
	float theta = 2.0*MX_PI*u;
	float phi=acos(2*v-1.0);
	float x=r*sin(phi)*sin(theta);
	float y=r*cos(phi)*sin(theta);
	float z=r*cos(theta);
	coordx.push_back(x);
	coordx.push_back(y);
	coordx.push_back(z);

	return coordx;
}

float randomAB(float A,float B)
{
	return A+(B-A)*rand()/(RAND_MAX);
 //   return A+(B-A)*randf0and1();
}

float Distance_point(vector<float> A,vector<float> B)
{
	float dis=0.0;
	for(int i=0;i<A.size();i++)
	{
		dis = dis + (A[i]-B[i])*(A[i]-B[i]);
	}
	dis=sqrt(dis);
	return dis;
}
//float Distance_pointy(vector<vector<float> > A,vector<vector<float> > B)
//{

//}
float Distance_pointx(vector<vector<float> > A,vector<vector<float> > B)
{
	float dis=0.0;
	float a1 = sqrt((A[0][0]-A[1][0])*(A[0][0]-A[1][0])+(A[0][1]-A[1][1])*(A[0][1]-A[1][1])+(A[0][2]-A[1][2])*(A[0][2]-A[1][2]));
	float a2 = sqrt((A[0][0]-A[2][0])*(A[0][0]-A[2][0])+(A[0][1]-A[2][1])*(A[0][1]-A[2][1])+(A[0][2]-A[2][2])*(A[0][2]-A[2][2]));
	float a3 = sqrt((A[0][0]-A[3][0])*(A[0][0]-A[3][0])+(A[0][1]-A[3][1])*(A[0][1]-A[3][1])+(A[0][2]-A[3][2])*(A[0][2]-A[3][2]));
	float b1 = sqrt((B[0][0]-B[1][0])*(B[0][0]-B[1][0])+(B[0][1]-B[1][1])*(B[0][1]-B[1][1])+(B[0][2]-B[1][2])*(B[0][2]-B[1][2]));
	float b2 = sqrt((B[0][0]-B[2][0])*(B[0][0]-B[2][0])+(B[0][1]-B[2][1])*(B[0][1]-B[2][1])+(B[0][2]-B[2][2])*(B[0][2]-B[2][2]));
	float b3 = sqrt((B[0][0]-B[3][0])*(B[0][0]-B[3][0])+(B[0][1]-B[3][1])*(B[0][1]-B[3][1])+(B[0][2]-B[3][2])*(B[0][2]-B[3][2]));

/*	for(int i=0;i<A.size();i++)
	{
		for(int j=0;j<A[i].size();j++)
		{
//			dis = dis + (A[i][j]-B[i][j])*(A[i][j]-B[i][j]);	
			dis = dis + (a1-b1)*(a1-b1)+		
		}
		
	} */
	dis = (a1-b1)*(a1-b1) + (a2-b2)*(a2-b2) + (a3-b3)*(a3-b3);
	dis=sqrt(dis);
	return dis;
} 

vector<float> calculateCoordinatesx(vector<float> refA,vector<float> refB,vector<float> refC,float L,float ang,float di)
{
	vector<float> AV(3,0.0),BV(3,0.0),CV(3,0.0);
	AV=refA;
	BV=refB;
	CV=refC;

	vector<float> CA(3,0.0),CB(3,0.0); 
	CA[0]=AV[0]-CV[0];
	CA[1]=AV[1]-CV[1];
	CA[2]=AV[2]-CV[2];
	CB[0]=BV[0]-CV[0];
	CB[1]=BV[1]-CV[1];
	CB[2]=BV[2]-CV[2];

	float AX = CA[0];
	float AY = CA[1];
	float AZ = CA[2];

	float BX = CB[0];
	float BY = CB[1];
	float BZ = CB[2];

	float A = (AY*BZ)-(AZ*BY);
	float B = (AZ*BX)-(AX*BZ);
	float G = (AX*BY)-(AY*BX);

	float F = sqrt(BX*BX + BY*BY + BZ*BZ) * L * cos(ang * (3.1415926/180.0));
//	float F = sqrt(BX*BX + BY*BY + BZ*BZ) * L * cos(ang);

	float constx = sqrt((pow(B*BZ-BY*G,2))*(-(F*F)*(A*A+B*B+G*G)+(B*B*(BX*BX+BZ*BZ) + A*A*(BY*BY+BZ*BZ)- (2*A*BX*BZ*G) + (BX*BX+ BY*BY)*G*G - (2*B*BY)*(A*BX+BZ*G))*L*L));
	float denom = (B*B)*(BX*BX+BZ*BZ)+ (A*A)*(BY*BY+BZ*BZ) - (2*A*BX*BZ*G) + (BX*BX+BY*BY)*(G*G) - (2*B*BY)*(A*BX+BZ*G);

	float X= ((B*B*BX*F)-(A*B*BY*F)+(F*G)*(-A*BZ+BX*G)+constx)/denom;

	float Y=0.0,Z=0.0;
    if((B==0 || BZ==0) && (BY==0 || G==0))
    {
        float const1=sqrt( G*G*(-A*A*X*X+(B*B+G*G)*(L-X)*(L+X)));
        Y= ((-A*B*X)+const1)/(B*B+G*G);
        Z= -(A*G*G*X+B*const1)/(G*(B*B+G*G));
    }
    else
    {
        Y= ((A*A*BY*F)*(B*BZ-BY*G)+ G*( -F*pow(B*BZ-BY*G,2) + BX*constx) - A*( B*B*BX*BZ*F- B*BX*BY*F*G + BZ*constx)) / ((B*BZ-BY*G)*denom);
        Z= ((A*A*BZ*F)*(B*BZ-BY*G) + (B*F)*pow(B*BZ-BY*G,2) + (A*BX*F*G)*(-B*BZ+BY*G) - B*BX*constx + A*BY*constx) / ((B*BZ-BY*G)*denom);
	}
	vector<float> D(3,0.0);
	D[0] = X + CV[0];
	D[1] = Y + CV[1];
	D[2] = Z + CV[2];

	float temp=0.0;

	temp=Points2Dihedral(AV,BV,CV,D) * (180.0/3.1415926);
//	cout<<"temp: "<<temp<<endl;

	di = di - temp;
	vector<float> Dx(3,0.0);
	Dx[0] = D[0] - BV[0];
	Dx[1] = D[1] - BV[1];
	Dx[2] = D[2] - BV[2];

	di = di * 3.1415926/180.0;
	GroupRotationt(CV,BV,di,Dx);
	D[0] = Dx[0] + BV[0];
	D[1] = Dx[1] + BV[1];
	D[2] = Dx[2] + BV[2];

	return D;
}

//struct distmap
//{
//    int ires;// i res
//    int jres; // j res
//    float dv; // dist value
//};

Model REMC_sampley(Model posex, string inputM,float MRC_R,float mapsp,string dist_name,char bindir[],char bindiry[])
{
    int nresobins = 200;  //  # resolution bins for statistics
    vector<float> resobins, mapI, maskedmapI, modelI, maskI;
    vector<int> resobin_counts;
    vector<float> perResCC, perResStrain;

    vector<float> mapmapFSC, maskedMapMapFSC;
    vector<float> modelmapFSC, maskedModelMapFSC;

    bool truncate_map = false ;
    bool maskonly = false ;
    bool cutonly = false ;
    bool bin_squared =  false;
    bool mask =false ;
    string alt_mapfile=" ";

    // resolution limits for analysis
    float hires = 0.0;  // high res limit
    float lowres = 1000.0;  // low res limit
    float truncate_hires = 0.0;  // high res truncation
    float truncate_lowres = 1000.0;  // low res truncation
    float mask_resolution = 0.0;   //  radius for masking
    bool perres =false;  //dump extra output

    ObjexxFCL::FArray3D< float > rhoC, rhoMask, rhoOmask, rhoO2mask;
    ObjexxFCL::FArray3D< std::complex<float> > FrhoC, FrhoMask, FrhoCmask, FrhoOmask, FrhoO2mask;
    ObjexxFCL::FArray3D< std::complex<float> > FrhoO, FrhoO2;   

    // read Density map
    ElectronDensity theDensityMap;
    theDensityMap.preso= MRC_R;
    theDensityMap.pATOM_MASK =7.0;
//    theDensityMap.pATOM_MASK =7.0;
    theDensityMap.pCA_MASK =7.0;
    theDensityMap.pforce_apix_on_map_load_= 0.0;
    theDensityMap.pWINDOW_ =1;
    theDensityMap.pscore_window_context_ = false ;
    theDensityMap.premap_symm_ = false ;
    theDensityMap.pde_edensity = false;
    theDensityMap.pnkbins_ = 0 ;   

    theDensityMap.readMRCandResize(inputM,MRC_R,mapsp);
    int Denx = theDensityMap.density.u1();
    int Deny = theDensityMap.density.u2();
    int Denz = theDensityMap.density.u3();
    int Denxyz = Denx*Deny*Denz;  
    vector<int> grid(3,0); 
    grid[0] = theDensityMap.grid[0];  
    grid[1] = theDensityMap.grid[1];
    grid[2] = theDensityMap.grid[2];
    vector<float> origin(3,0.0);
    origin[0] = theDensityMap.origin[0]; 
    origin[1] = theDensityMap.origin[1];
    origin[2] = theDensityMap.origin[2];    
    float padding = theDensityMap.ATOM_MASK_PADDING; 
    float ca_m = theDensityMap.CA_MASK;
    float atom_m = theDensityMap.ATOM_MASK; 
//  theDensityMap.calcRhoC(posey,4.0,rhoC,rhoMask);
//    theDensityMap.calcRhoCx(posex,4.0,rhoC,rhoMask);
//  theDensityMap.calcRhoCy(pose,4.0,rhoC,rhoMask);
//    float CC;
//  CC=theDensityMap.getRSCC(rhoC,rhoMask);
 //   CC=theDensityMap.getRSCCX(rhoC,rhoMask);
 //   cout<<"CC: "<<CC<<endl; 
    cout<<"origin: "<<theDensityMap.origin[0]<<" "<<theDensityMap.origin[1]<<" "<<theDensityMap.origin[2]<<endl;;

    clock_t start,finish;
    double totaltime;
    start=clock(); 

    Model fin_pose;
//    fin_pose = posex;
//  vector<Rot> fin_vec;
//  vector<Rot> beg_vec;
//  float min_dis=100.0;
//  vector<vector<float> > point;
//  vector<float> old_v(3,0.0),old_v1(3,0.0);
//  vector<vector<float> > old_v,old_v1,old_tmp,new_v;
//  old_v=vector<vector<float> >(4,vector<float>(3,0.0));
//  old_v1=vector<vector<float> >(4,vector<float>(3,0.0));
//  old_tmp = vector<vector<float> >(4,vector<float>(3,0.0));
//  new_v = vector<vector<float> >(4,vector<float>(3,0.0));
//  vector<float> new_v(3,0.0);
    float KT=0.001;
//  vector<vector<Rot> > frag_pre(7),frag_nat(7);


    
    vector<Rot> pointsx;
//  vector<vector<float> > pointsABx; // all backbone atom
    vector<poseCoord> pointsBx;
//  vector<vector<float> > pointsABx0; 
    vector<poseCoord> pointsBx0;
    vector<string> pointsrenm;
//  vector<vector<float> > pointsx,pointsy;
//  int pnum = posex.size();

//  cout<<"33 "<<endl;
//  cout<< posex.size()<<endl;
//  cout<< posey.size()<<endl;
//  for(int j=0;j<pnum;j++)
//  {
//      pointsx.push_back(posex[j].x_);
//      pointsy.push_back(posey[j].x_);
//  }

    int pnum =0,pnumx=0;
//  cout<< " 1111"<<endl;
    for(int i=0;i<posex.chains.size();i++)
    {
//      cout<< "chainx size: "<<posex.chains.size()<<endl;
//      cout<< "chainy size: "<<posey.chains.size()<<endl;
        Chain Chanx = posex.chains[i];
//      Chain Chany = posey.chains[i];
        for(int j=0;j<Chanx.residues.size();j++)
        {
//          cout<< "residuex size: "<< Chanx.residues.size()<<endl;
//          cout<< "residuey size: "<< Chany.residues.size()<<endl;
            Residue Resdx = Chanx.residues[j];
    //      Residue Resdy = Chany.residues[j];
            Rot rotx;
//          cout<<"Res: "<<Resdx.atoms.size()<<endl;
            pointsrenm.push_back(Resdx.resname);
//          pointsrenmy.push_back(Resdy.resname);
            for(int t=0;t<Resdx.atoms.size();t++)
            {
//              cout<< "atomx size: "<<Resdx.atoms.size()<<endl;
                poseCoord pcdx;
                Atom Atmx=Resdx.atoms[t];
                if(Atmx.atona != " CA " && Atmx.atona != " C  " && Atmx.atona != " N  " && Atmx.atona != " O  ")
                {
                    rotx.res_atm.push_back(Atmx.xyzVect);
                    vector<string> atmnmx = string_splitxx(Atmx.atona,' ');
                    string atmnmy = atmnmx[0];
                    rotx.resatmnm.push_back(atmnmy);
                }               
                if(Atmx.atona == " CA ") {rotx.ca = Atmx.xyzVect;pcdx.x_=Atmx.xyzVect;pcdx.B_=Atmx.tempfac;pcdx.elt_= name2eltx(Atmx.atona);pointsBx.push_back(pcdx);pnumx = pnumx +1;}
                if(Atmx.atona == " C  ") {rotx.cc = Atmx.xyzVect;pcdx.x_=Atmx.xyzVect;pcdx.B_=Atmx.tempfac;pcdx.elt_= name2eltx(Atmx.atona);pointsBx.push_back(pcdx);}//pnumx = pnumx +1;
                if(Atmx.atona == " N  ") {rotx.cn = Atmx.xyzVect;pcdx.x_=Atmx.xyzVect;pcdx.B_=Atmx.tempfac;pcdx.elt_= name2eltx(Atmx.atona);pointsBx.push_back(pcdx);}//pnumx = pnumx +1;
            /*    if(Atmx.atona == " C  ") {rotx.cc = Atmx.xyzVect;}
                if(Atmx.atona == " N  ") {rotx.cn = Atmx.xyzVect;}
                if(Atmx.atona == " O  ") {rotx.co = Atmx.xyzVect;}
                if(Atmx.atona == " CB ") {rotx.cb = Atmx.xyzVect;}
                if(Atmx.atona == " CG ") {rotx.cg = Atmx.xyzVect;}
                if(Atmx.atona == " CG1") {rotx.cg1 = Atmx.xyzVect;}
                if(Atmx.atona == " CG2") {rotx.cg2 = Atmx.xyzVect;}
                if(Atmx.atona == " OG ") {rotx.og  = Atmx.xyzVect;}
                if(Atmx.atona == " SG ") {rotx.sg = Atmx.xyzVect;}
                if(Atmx.atona == " CD ") {rotx.cd = Atmx.xyzVect;}
                if(Atmx.atona == " SD ") {rotx.sd = Atmx.xyzVect;}
                if(Atmx.atona == " CD1") {rotx.cd1 = Atmx.xyzVect;}
                if(Atmx.atona == " CD2") {rotx.cd2 = Atmx.xyzVect;}
                if(Atmx.atona == " ND2") {rotx.nd2 = Atmx.xyzVect;}
                if(Atmx.atona == " ND1") {rotx.nd1 = Atmx.xyzVect;}
                if(Atmx.atona == " OD2") {rotx.od2 = Atmx.xyzVect;}
                if(Atmx.atona == " OE2") {rotx.oe2 = Atmx.xyzVect;}
                if(Atmx.atona == " CE ") {rotx.ce = Atmx.xyzVect;}
                if(Atmx.atona == " NE2") {rotx.ne2 = Atmx.xyzVect;}
                if(Atmx.atona == " NE ") {rotx.ne = Atmx.xyzVect;}
                if(Atmx.atona == " NZ ") {rotx.nz = Atmx.xyzVect;}
                if(Atmx.atona == " CZ ") {rotx.cz = Atmx.xyzVect;}             */     
    /*          if(Atmy.atona == " CB ") roty.cb = Atmy.xyzVect;                    
                if(Atmx.atona == " CG ") rotx.cg = Atmy.xyzVect;
                if(Atmy.atona == " CG ") roty.cg = Atmy.xyzVect;
                if(Atmx.atona == " CD ") rotx.cd = Atmy.xyzVect;
                if(Atmy.atona == " CD ") roty.cd = Atmy.xyzVect;
                if(Atmx.atona == " CE ") rotx.ce = Atmy.xyzVect;
                if(Atmy.atona == " CE ") roty.ce = Atmy.xyzVect;
                if(Atmx.atona == " CZ ") rotx.cz = Atmy.xyzVect;
                if(Atmy.atona == " CZ ") roty.cz = Atmy.xyzVect;        */                  

            }
    /*      for(int t=0;t<Resdy.atoms.size();t++)
            {
//              cout<< "atomy size: "<<Resdy.atoms.size()<<endl;
                Atom Atmy=Resdy.atoms[t];
                if(Atmy.atona != " CA " && Atmy.atona != " C  " && Atmy.atona != " N  " && Atmy.atona != " O  ")
                {
                    roty.res_atm.push_back(Atmy.xyzVect);
                }               
                if(Atmy.atona == " CA ") {roty.ca = Atmy.xyzVect;pointsABy.push_back(Atmy.xyzVect);}
                if(Atmy.atona == " C  ") {roty.cc = Atmy.xyzVect;pointsABy.push_back(Atmy.xyzVect);}
                if(Atmy.atona == " N  ") {roty.cn = Atmy.xyzVect;pointsABy.push_back(Atmy.xyzVect);}
                if(Atmy.atona == " O  ") {roty.co = Atmy.xyzVect;}
                if(Atmy.atona == " CB ") {roty.cb = Atmy.xyzVect;}
                if(Atmy.atona == " CG ") {roty.cg = Atmy.xyzVect;}
                if(Atmy.atona == " CG1") {roty.cg1 = Atmy.xyzVect;}
                if(Atmy.atona == " CG2") {roty.cg2 = Atmy.xyzVect;}
                if(Atmy.atona == " OG ") {roty.og  = Atmy.xyzVect;}
                if(Atmy.atona == " SG ") {roty.sg = Atmy.xyzVect;}
                if(Atmy.atona == " CD ") {roty.cd = Atmy.xyzVect;}
                if(Atmy.atona == " SD ") {roty.sd = Atmy.xyzVect;}
                if(Atmy.atona == " CD1") {roty.cd1 = Atmy.xyzVect;}
                if(Atmy.atona == " CD2") {roty.cd2 = Atmy.xyzVect;}
                if(Atmy.atona == " ND2") {roty.nd2 = Atmy.xyzVect;}
                if(Atmy.atona == " ND1") {roty.nd1 = Atmy.xyzVect;}
                if(Atmy.atona == " OD2") {roty.od2 = Atmy.xyzVect;}
                if(Atmy.atona == " OE2") {roty.oe2 = Atmy.xyzVect;}
                if(Atmy.atona == " CE ") {roty.ce = Atmy.xyzVect;}
                if(Atmy.atona == " NE2") {roty.ne2 = Atmy.xyzVect;}
                if(Atmy.atona == " NE ") {roty.ne = Atmy.xyzVect;}
                if(Atmy.atona == " NZ ") {roty.nz = Atmy.xyzVect;}
                if(Atmy.atona == " CZ ") {roty.cz = Atmy.xyzVect;}              
            }           */
            pointsx.push_back(rotx);
//          pointsy.push_back(roty);
            pnum = pnum +1; 
//          cout<< pnum<<endl;      
        }
    }
    cout<<"pnum,pnumx: "<<pnum<<" "<<pnumx<<endl;
    pointsBx0 = pointsBx;

    cout<<"binddie: "<<bindir<<endl;
    char oneline[600];
    sprintf(oneline,"%s/bnramaprob",bindir);
    vector<vector<vector<float>>> rlogduke;
    vector<vector<vector<int>>> ramaduke;    
    loadramadukebn(oneline,rlogduke,ramaduke); // torsion angle distribution
//    loadramadukebn(oneline);    
//    sprintf(oneline,"%s/newsgdistriaanew72.txt",bindir); // side-chane 
//    vector<vector<vector<vector<float>>>> sgpos2;    
//    loadsgpos2(oneline,72,sgpos2);    


//    cout<<"SS"<<endl;
/*
    vector<distmap> dist_map;
    for(int i=0;i<pnum;i++)
    {
        for(int j=i;j<pnum;j++)
        {
            if(abs(i-j)<2) continue;
            float tmp_dist = Distance_point(pointsBx[3*i+1].x_,pointsBx[3*j+1].x_);
            if(tmp_dist<6.5)
            {
                distmap tmp_map;
                tmp_map.ires = i;
                tmp_map.jres = j;
                tmp_map.dv = tmp_dist;
                dist_map.push_back(tmp_map);
            }
        }
    } */
//    int N_dist=dist_pre.size();
    vector<float> axyz;
    axyz=vector<float>(3,0.0);
    float ang; 
    float angle_rotate; 
    float old_dE1=0.0;
    float new_dE1=0.0;    
    float old_CC=0.0;
    float new_CC=0.0;
    float old_dE2=0.0;
    float new_dE2=0.0; 
    float old_vwd1=0.0;
    float old_vwd2=0.0;
    float new_vwd1=0.0;
    float new_vwd2=0.0; 
    float new_dE=0.0;
    float old_dE=0.0;       
    float dE=0.0;
    float old_vwd= 0.0;
    float new_vwd= 0.0;
    float tmp_den=0.0;
    float tmp_rhc=0.0;
    float old_clash = 0.0;
    float new_clash = 0.0;
//    vector<poseCoord> fin_mat;
//    vector<vector<float> > fin_matx;
    vector<float> coord_change(3,0.0);

    vector<float> atm_idx(3,0.0);


    float effReso = std::max( 2.4+0.8*MRC_R , double(MRC_R ));
    float k=(M_PI/effReso)*(M_PI/effReso);
    float a=33.0;  // treat everything as ALA
    float C=a*pow(k/3.1415926,1.5);
    int tmp_i=0;
    ObjexxFCL::FArray3D< float > rhoC0;
//    ObjexxFCL::FArray3D< float > rhoC01;
    ObjexxFCL::FArray3D< float > inv_rho_mask0;
//    ObjexxFCL::FArray3D< float > inv_rho_mask1;
    rhoC0.dimension(theDensityMap.density.u1() , theDensityMap.density.u2() , theDensityMap.density.u3());
//    rhoC01.dimension(theDensityMap.density.u1() , theDensityMap.density.u2() , theDensityMap.density.u3());
    inv_rho_mask0.dimension(theDensityMap.density.u1() , theDensityMap.density.u2() , theDensityMap.density.u3());
//    inv_rho_mask1.dimension(theDensityMap.density.u1() , theDensityMap.density.u2() , theDensityMap.density.u3());
    for ( int t=0; t<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++t ) {
        rhoC0[t]=0.0;
        inv_rho_mask0[t]=1.0;
    }
    vector<float> del_ijx(3,0.0);
    vector<float> atm_jx(3,0.0); 
    string elt_i; 
    float bond0max=0.0,bond1max=0.0,bond2max=0.0;
    float bond0min=1000.0,bond1min=1000.0,bond2min=1000.0;        
    for(int i=0; i<pointsBx.size();i++)
    {
//        cout<<"YYYY: "<<pointsBx[i].elt_<<endl;;
//        if(pointsBx[i].elt_ == "CA")
//        {
         if(pointsBx[i].elt_ == "CA")
        {           
//            cout<<"RRRR: "<<i;
            vector<float> cartX1;
            vector<float> fracX1;
            elt_i = pointsBx[i].elt_;
            elt_i = elt_i[0];
            OneGaussianScattering sig_j = get_A( elt_i );
            k = sig_j.k( theDensityMap.effectiveB );
            C = sig_j.C( k );
            if ( C < 1e-6 ) continue;   

            cartX1 = pointsBx[i].x_;           
            MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);

            // the location of atom in grid ?
            atm_idx[0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
            atm_idx[1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
            atm_idx[2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);   
    //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
    //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
    //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
    //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
            for(int z=1;z<=theDensityMap.density.u3();z++)
            {
                atm_jx[2] = z;
                del_ijx[2] =(atm_idx[2]-atm_jx[2])/theDensityMap.grid[2];
                if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                del_ijx[0] = del_ijx[1] = 0.0;
                vector<float> frac_tmpz;
                MatrixTimesTransVector(theDensityMap.f2c,del_ijx,frac_tmpz);
                if(square_len(frac_tmpz)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                    
                for(int y=1;y<=theDensityMap.density.u2();y++)
                {
                    atm_jx[1] = y;
                    del_ijx[1] = (atm_idx[1] - atm_jx[1])/theDensityMap.grid[1];
                    // wrap-around??
                    if ( del_ijx[1] > 0.5 ) del_ijx[1]-=1.0;
                    if ( del_ijx[1] < -0.5 ) del_ijx[1]+=1.0;                        
                    del_ijx[0] = 0.0;
                    vector<float> frac_tmpy;
                    MatrixTimesTransVector(theDensityMap.f2c,del_ijx,frac_tmpy);
                    if(square_len(frac_tmpy)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                        
                    for(int x=1;x<=theDensityMap.density.u1();x++)
                    {
                        atm_jx[0] = x;
                        del_ijx[0] = (atm_idx[0] - atm_jx[0])/theDensityMap.grid[0];
                        // wrap-around??
                        if ( del_ijx[0] > 0.5 ) del_ijx[0]-=1.0;
                        if ( del_ijx[0] < -0.5 ) del_ijx[0]+=1.0;                            
                        vector<float> cart_del_ij2;
                        MatrixTimesTransVector(theDensityMap.f2c,del_ijx,cart_del_ij2);
                        float d2 = square_len(cart_del_ij2);
                        if(d2 > (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;
                    
                        float atm = C*exp(-k*d2);
//                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
                        float sigmoid_msk = exp( d2 - (theDensityMap.ATOM_MASK)*(theDensityMap.ATOM_MASK)  );
                    //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
                    //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
                        float inv_msk = 1/(1+sigmoid_msk);
                        rhoC0(x,y,z) += atm;
                        inv_rho_mask0(x,y,z) *= (1 - inv_msk);                       

//                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
//                                mask2(x,y,z) = 1.0; // problem?
//                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
//                                    float atm = C*exp(-k*d2);
//                                    rhoC(x,y,z) += atm;
//                                }
//                            }
                    }
                }
            } 
//            cout<<" "<<tmp_i<<" "<<endl;           
            tmp_i = tmp_i + 1;
        }    
//        }
    }   
    float max_den = 0.0;
    float max_rhc = 0.0;
    for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
//        cout<<theDensityMap.density[x]<<" ";
        if(theDensityMap.density[x]>max_den) max_den = theDensityMap.density[x] ;
        if(rhoC0[x]>max_rhc) max_rhc =rhoC0[x] ;
    } 
    vector<float> mass_pdb(3,0.0);
    vector<float> mass_mrc(3,0.0);
    vector<float> coor_pdb(3,0.0);  
    float allpdbmass=0.0;
    float allmrcmass=0.0;
    float sumC_i2=0.0, sumO_i2=0.0, sumCO_i2=0.0, vol_i2=0.0, CC_i2=0.0;
    float sumO2_i2=0.0, sumC2_i2=0.0, varC_i2=0.0, varO_i2=0.0;
    float clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;    
    mass_pdb=vector<float>(3,0.0);
    mass_mrc=vector<float>(3,0.0);            
    for(int x=1;x<=theDensityMap.density.u1();x++)
    {
       for(int y=1;y<=theDensityMap.density.u2();y++)
        {
            for(int z=1;z<=theDensityMap.density.u3();z++)
            {
                clc_x2 = rhoC0(x,y,z);
//                clc_x2 = tmp_rhc;
                obs_x2 = theDensityMap.density(x,y,z);
//                obs_x2 = tmp_den;
                eps_x2 = 1.0-inv_rho_mask0(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
//                sumCO_i2 += eps_x2*clc_x2*obs_x2;
//                sumO_i2  += eps_x2*obs_x2;
//                sumO2_i2 += eps_x2*obs_x2*obs_x2;
//                sumC_i2  += eps_x2*clc_x2;
//                sumC2_i2 += eps_x2*clc_x2*clc_x2;
//                vol_i2   += eps_x2;
//                mass_pdb[0]=  mass_pdb[0] + eps_x2*clc_x2*x;
//                mass_pdb[1]=  mass_pdb[1] + eps_x2*clc_x2*y;
//                mass_pdb[2]=  mass_pdb[2] + eps_x2*clc_x2*z;
                mass_mrc[0]=  mass_mrc[0] + eps_x2*obs_x2*x;
                mass_mrc[1]=  mass_mrc[1] + eps_x2*obs_x2*y;
                mass_mrc[2]=  mass_mrc[2] + eps_x2*obs_x2*z;
                allpdbmass = allpdbmass + eps_x2*clc_x2;
                allmrcmass = allmrcmass + eps_x2*obs_x2;          
/*                mass_pdb[0]=  mass_pdb[0] + clc_x2*x;
                mass_pdb[1]=  mass_pdb[1] + clc_x2*y;
                mass_pdb[2]=  mass_pdb[2] + clc_x2*z;
                mass_mrc[0]=  mass_mrc[0] + obs_x2*x;
                mass_mrc[1]=  mass_mrc[1] + obs_x2*y;
                mass_mrc[2]=  mass_mrc[2] + obs_x2*z;
                allpdbmass = allpdbmass + clc_x2;
                allmrcmass = allmrcmass + obs_x2;   */                            
            }
        }
    } 
//    mass_pdb[0] = mass_pdb[0]/allpdbmass;
//    mass_pdb[1] = mass_pdb[1]/allpdbmass;
//    mass_pdb[2] = mass_pdb[2]/allpdbmass;
    mass_mrc[0] = mass_mrc[0]/allmrcmass;
    mass_mrc[1] = mass_mrc[1]/allmrcmass;
    mass_mrc[2] = mass_mrc[2]/allmrcmass; 
    cout<<"mass_mrc: "<<mass_mrc[0]<<" "<<mass_mrc[1]<<" "<<mass_mrc[2]<<endl;
    mass_mrc[0] = theDensityMap.density.u1()/2.0+(theDensityMap.origin[0]+1.0);//*(1.0-1.0/2.0);
    mass_mrc[1] = theDensityMap.density.u2()/2.0+(theDensityMap.origin[1]+1.0);//*(1.0-1.0/2.0);
    mass_mrc[2] = theDensityMap.density.u3()/2.0+(theDensityMap.origin[2]+1.0);//*(1.0-1.0/2.0);    
    cout<<"mass_mrc: "<<mass_mrc[0]<<" "<<mass_mrc[1]<<" "<<mass_mrc[2]<<endl;     
//    float vol_s = theDensityMap.voxel_volume();
//    vol_s = pow(vol_s,1.0/3.0);
//    mass_mrc[0]=mass_mrc[0]*vol_s;
//    mass_mrc[1]=mass_mrc[1]*vol_s;   
//    mass_mrc[2]=mass_mrc[2]*vol_s;
//    mass_mrc[0]=mass_mrc[0]*vol_s*theDensityMap.grid[0]/(theDensityMap.grid[0]+theDensityMap.grid[1]+theDensityMap.grid[2]);    
//    mass_mrc[1]=mass_mrc[1]*vol_s*theDensityMap.grid[1]/(theDensityMap.grid[0]+theDensityMap.grid[1]+theDensityMap.grid[2]);     
//    mass_mrc[2]=mass_mrc[2]*vol_s*theDensityMap.grid[2]/(theDensityMap.grid[0]+theDensityMap.grid[1]+theDensityMap.grid[2]);    
//    cout<<"vol_s: "<<vol_s<<endl;   
    mass_mrc[0] = mass_mrc[0]/theDensityMap.grid[0];
    mass_mrc[1] = mass_mrc[1]/theDensityMap.grid[1];
    mass_mrc[2] = mass_mrc[2]/theDensityMap.grid[2];
//    float tran0=Distance_point(mass_pdb,mass_mrc); 
//    vector<float> tran(3,0.0);
//    tran[0]=mass_pdb[0]- mass_mrc[0];
//    tran[1]=mass_pdb[1]- mass_mrc[1];
//    tran[2]=mass_pdb[2]- mass_mrc[2];
//    tran[0]=(tran[1] + (theDensityMap.origin[0]+1)*(1.0))/theDensityMap.grid[0];
//    tran[1]=(tran[1] + (theDensityMap.origin[0]+1)*(1.0))/theDensityMap.grid[0];
//    tran[2]=(tran[1] + (theDensityMap.origin[0]+1)*(1.0))/theDensityMap.grid[0];
//    cout<<"tran: "<<tran[0]<<" "<<tran[1]<<" "<<tran[2]<<endl;
    vector<float> mass_mrc0(3,0.0);
    MatrixTimesTransVector(theDensityMap.f2c,mass_mrc,mass_mrc0);
//    mass_mrc0[0] = mass_mrc0[0] + theDensityMap.origin[0];
//    mass_mrc0[1] = mass_mrc0[1] + theDensityMap.origin[1];
//    mass_mrc0[2] = mass_mrc0[2] + theDensityMap.origin[2];    
//    mass_mrc0[0]=mass_mrc0[0]/theDensityMap.grid[0];
//    mass_mrc0[1]=mass_mrc0[1]/theDensityMap.grid[1];
//    mass_mrc0[2]=mass_mrc0[2]/theDensityMap.grid[2];
//    mass_mrc0[0] = mass_mrc0[0]+ (theDensityMap.origin[0]+1.0);
//    mass_mrc0[1] = mass_mrc0[1]+ (theDensityMap.origin[1]+1.0);
//    mass_mrc0[2] = mass_mrc0[2]+ (theDensityMap.origin[2]+1.0);
    cout<<"mass_mrc0: "<<mass_mrc0[0]<<" "<<mass_mrc0[1]<<" "<<mass_mrc0[2]<<endl;
    coor_pdb=vector<float>(3,0.0);   
    for(int j=0;j<pnum;j++)
    {        
        coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
        coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
        coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];       
    }
    coor_pdb[0] = coor_pdb[0]/(pnum);
    coor_pdb[1] = coor_pdb[1]/(pnum);
    coor_pdb[2] = coor_pdb[2]/(pnum);
    cout<<"coor_pdb: "<<coor_pdb[0]<<" "<<coor_pdb[1]<<" "<<coor_pdb[2]<<endl;        
    vector<float> tran0(3,0.0);
    tran0[0]= mass_mrc0[0] - coor_pdb[0];
    tran0[1]= mass_mrc0[1] - coor_pdb[1];
    tran0[2]= mass_mrc0[2] - coor_pdb[2];    
    cout<<"tran0: "<<tran0[0]<<" "<<tran0[1]<<" "<<tran0[2]<<endl;
    for(int j=0;j<pointsBx.size();j++) // not last N
    {
//        pointsBx[j].x_[0] = tran0[0] + pointsBx[j].x_[0];
//        pointsBx[j].x_[1] = tran0[1] + pointsBx[j].x_[1];
//        pointsBx[j].x_[2] = tran0[2] + pointsBx[j].x_[2];
//        tmp3=tmp3+1;
    }  

    ElectronDensity theDensityMapx;
    theDensityMapx.preso= MRC_R;
    theDensityMapx.pATOM_MASK =7.0;
    theDensityMapx.pCA_MASK =7.0;
    theDensityMapx.pforce_apix_on_map_load_= 0.0;
    theDensityMapx.pWINDOW_ =1;
    theDensityMapx.pscore_window_context_ = false ;
    theDensityMapx.premap_symm_ = false ;
    theDensityMapx.pde_edensity = false;
    theDensityMapx.pnkbins_ = 0 ;  
    float mapspx=2.0;
    float new_CC1=0.0;
    float old_CC1=0.0;
    float new_CC2=0.0;     
    theDensityMapx.readMRCandResize(inputM,MRC_R,mapspx);
    int pDenx = theDensityMapx.density.u1();
    int pDeny = theDensityMapx.density.u2();
    int pDenz = theDensityMapx.density.u3();
    int pDenxyz = pDenx*pDeny*pDenz;  
    vector<int> pgrid(3,0); 
    pgrid[0] = theDensityMapx.grid[0];  
    pgrid[1] = theDensityMapx.grid[1];
    pgrid[2] = theDensityMapx.grid[2];
    vector<float> porigin(3,0.0);
    porigin[0] = theDensityMapx.origin[0]; 
    porigin[1] = theDensityMapx.origin[1];
    porigin[2] = theDensityMapx.origin[2];    
    float ppadding = theDensityMapx.ATOM_MASK_PADDING; 
    float pca_m = theDensityMapx.CA_MASK;
    float patom_m = theDensityMapx.ATOM_MASK;

    vector<poseCoord> best_model;
    float best_E=10.0;
    ObjexxFCL::FArray3D< float > rhoC0x;
    ObjexxFCL::FArray3D< float > inv_rho_mask0x;
    rhoC0x.dimension(theDensityMapx.density.u1() , theDensityMapx.density.u2() , theDensityMapx.density.u3());
    inv_rho_mask0x.dimension(theDensityMapx.density.u1() , theDensityMapx.density.u2() , theDensityMapx.density.u3());
    // find initial position for A and B
    float init0max=-10000.0,init1max=-10000.0,init2max=-10000.0;
    float init0min=10000.0,init1min=10000.0,init2min=10000.0;       
    best_E=10.0;
    coor_pdb=vector<float>(3,0.0);   
    for(int j=0;j<pnum;j++)
    {        
        coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
        coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
        coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];       
    }
    coor_pdb[0] = coor_pdb[0]/(pnum);
    coor_pdb[1] = coor_pdb[1]/(pnum);
    coor_pdb[2] = coor_pdb[2]/(pnum);      
    for ( int t=0; t<pDenxyz; ++t ) {
        rhoC0x[t]=0.0;
        inv_rho_mask0x[t]=1.0;
    }  
    del_ijx=vector<float>(3,0.0);
    atm_jx=vector<float>(3,0.0);
    tmp_i=0; 
//    string elt_i; 
//    float bond0max=0.0,bond1max=0.0,bond2max=0.0;
//    float bond0min=1000.0,bond1min=1000.0,bond2min=1000.0;   
    bond0max=-10000.0;bond1max=-10000.0;bond2max=-10000.0;
    bond0min=10000.0;bond1min=10000.0;bond2min=10000.0;       
    for(int i=0; i<pointsBx.size();i++)
    {
//        cout<<"YYYY: "<<pointsBx[i].elt_<<endl;;
//        if(pointsBx[i].elt_ == "CA")
//        {
        if(pointsBx[i].elt_ == "CA")
        {           
//            cout<<"RRRR: "<<i;
            vector<float> cartX1;
            vector<float> fracX1;
            elt_i = pointsBx[i].elt_;
            elt_i = elt_i[0];
            OneGaussianScattering sig_j = get_A( elt_i );
            k = sig_j.k( theDensityMapx.effectiveB );
            C = sig_j.C( k );
            if ( C < 1e-6 ) continue;   

            cartX1 = pointsBx[i].x_;           
            MatrixTimesTransVector(theDensityMapx.c2f,cartX1,fracX1);

            // the location of atom in grid ?
            atm_idx[0] = pos_mod (double(fracX1[0]*theDensityMapx.grid[0] - theDensityMapx.origin[0] + 1) , (double)theDensityMapx.grid[0]);
            atm_idx[1] = pos_mod (double(fracX1[1]*theDensityMapx.grid[1] - theDensityMapx.origin[1] + 1) , (double)theDensityMapx.grid[1]);
            atm_idx[2] = pos_mod (double(fracX1[2]*theDensityMapx.grid[2] - theDensityMapx.origin[2] + 1) , (double)theDensityMapx.grid[2]);   
    //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
    //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
    //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
    //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
            for(int z=1;z<=theDensityMapx.density.u3();z++)
            {
                atm_jx[2] = z;
                del_ijx[2] =(atm_idx[2]-atm_jx[2])/theDensityMapx.grid[2];
                if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                del_ijx[0] = del_ijx[1] = 0.0;
                vector<float> frac_tmpz;
                MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,frac_tmpz);
                if(square_len(frac_tmpz)> (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;                    
                for(int y=1;y<=theDensityMapx.density.u2();y++)
                {
                    atm_jx[1] = y;
                    del_ijx[1] = (atm_idx[1] - atm_jx[1])/theDensityMapx.grid[1];
                    // wrap-around??
                    if ( del_ijx[1] > 0.5 ) del_ijx[1]-=1.0;
                    if ( del_ijx[1] < -0.5 ) del_ijx[1]+=1.0;                        
                    del_ijx[0] = 0.0;
                    vector<float> frac_tmpy;
                    MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,frac_tmpy);
                    if(square_len(frac_tmpy)> (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;                        
                    for(int x=1;x<=theDensityMapx.density.u1();x++)
                    {
                        atm_jx[0] = x;
                        del_ijx[0] = (atm_idx[0] - atm_jx[0])/theDensityMapx.grid[0];
                        // wrap-around??
                        if ( del_ijx[0] > 0.5 ) del_ijx[0]-=1.0;
                        if ( del_ijx[0] < -0.5 ) del_ijx[0]+=1.0;                            
                        vector<float> cart_del_ij2;
                        MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,cart_del_ij2);
                        float d2 = square_len(cart_del_ij2);
                        if(d2 > (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;
                    
                        float atm = C*exp(-k*d2);
//                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
                        float sigmoid_msk = exp( d2 - (theDensityMapx.ATOM_MASK)*(theDensityMapx.ATOM_MASK)  );
                    //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
                    //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
                        float inv_msk = 1/(1+sigmoid_msk);
                        rhoC0x(x,y,z) += atm;
                        inv_rho_mask0x(x,y,z) *= (1 - inv_msk);

                        if(x>bond0max) bond0max = x;
                        if(y>bond1max) bond1max = y;
                        if(z>bond2max) bond2max = z;
                        if(x<bond0min) bond0min = x;
                        if(y<bond1min) bond1min = y;
                        if(z<bond2min) bond2min = z;                                               

//                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
//                                mask2(x,y,z) = 1.0; // problem?
//                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
//                                    float atm = C*exp(-k*d2);
//                                    rhoC(x,y,z) += atm;
//                                }
//                            }
                    }
                }
            } 
//            cout<<" "<<tmp_i<<" "<<endl;           
            tmp_i = tmp_i + 1;
        }    
//        }
    }        
    sumC_i2=0.0; sumO_i2=0.0; sumCO_i2=0.0; vol_i2=0.0; CC_i2=0.0;
    sumO2_i2=0.0; sumC2_i2=0.0; varC_i2=0.0; varO_i2=0.0;
    clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
//    allpdbmass=0.0;
//    allmrcmass=0.0;
//    mass_pdb=vector<float>(3,0.0);
//    mass_mrc=vector<float>(3,0.0);
    for(int x=int(bond0min);x<=bond0max;x++)
    {
       for(int y=int(bond1min);y<=bond1max;y++)
        {
            for(int z=int(bond2min);z<=bond2max;z++)
            {
                clc_x2 = rhoC0x(x,y,z);
//                clc_x2 = tmp_rhc;
                obs_x2 = theDensityMapx.density(x,y,z);
//                obs_x2 = tmp_den;
                eps_x2 = 1-inv_rho_mask0x(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                sumCO_i2 += eps_x2*clc_x2*obs_x2;
                sumO_i2  += eps_x2*obs_x2;
                sumO2_i2 += eps_x2*obs_x2*obs_x2;
                sumC_i2  += eps_x2*clc_x2;
                sumC2_i2 += eps_x2*clc_x2*clc_x2;
                vol_i2   += eps_x2;
            }
        }
    } 
/*   for ( int x=0; x<theDensityMapx.density.u1()*theDensityMapx.density.u2()*theDensityMapx.density.u3(); ++x ) {
        clc_x2 = rhoC0x[x];
//            clc_x2 = tmp_rhc;
        obs_x2 = theDensityMapx.density[x];
//            obs_x2 = tmp_den;
        eps_x2 = 1-inv_rho_mask0x[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
//            eps_x2 = 1.0;

        // SMOOTHED
        sumCO_i2 += eps_x2*clc_x2*obs_x2;
        sumO_i2  += eps_x2*obs_x2;
        sumO2_i2 += eps_x2*obs_x2*obs_x2;
        sumC_i2  += eps_x2*clc_x2;
        sumC2_i2 += eps_x2*clc_x2*clc_x2;
        vol_i2   += eps_x2;
    }    */
    varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
    varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 );
//    cout<<"sumO2_i2: "<<sumO2_i2<<" "<<sumO_i2<<" "<<vol_i2<<endl;
//    cout<<"varC_i2: "<<varC_i2<<" "<<varO_i2<<endl;
    if ( varC_i2 == 0 || varO_i2 == 0 ) {
        CC_i2 = 0;
    } else {
        CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
    }   
    new_CC2=1.0-CC_i2;
//    cout<<"new_CC2: "<<new_CC2<<endl;
    for(int iii=0;iii<0;iii++) // 500
    {
//        cout<<"trdisx,anglex,KT: "<<trdis0n<<" "<<angle0n<<" "<<KT<<endl;
//        float asin_theta=2.0*RandomDoubleX(0.0,1.0)-1.0;
        float asin_theta=2.0*randf0and1()-1.0;
        float acos_theta=sqrt(1.0-asin_theta*asin_theta);
        float apha=2.0*PI*randf0and1();
        float awx=acos_theta*cos(apha);
        float awy=acos_theta*sin(apha);
        float awz=asin_theta;
        // Translation Vector
        float t0=1.0;
    //    float t0=1.0;

    //    cout<<"t0: "<<t1<<" "<<t2<<" "<<t3<<endl;

        // Rotation matrix
        float anggg=60.0;
        float angle_rotategg=(2.0*randf0and1()-1.0)*anggg; // rotate angle
    //    if(new_CC2>0.95) t0=3.0;
    /*    if(acc_rate>40 && acc_rate<50)
        {
//            acc_rate=0;
        //    anggg=float(rand()%90)+120.0;
            anggg=float(randIntCustom(0,90))+90.0;
       //     anggg=float(randIntCustom(0,120));
            float tmp_rv=(2.0*randf0and1()-1.0);
            if(tmp_rv>=0.0)
            {
              angle_rotategg=anggg; // rotate angle  
            }
            else
            {
                angle_rotategg=(-1.0)*anggg; // rotate angle  
            }
        //    t0=float(rand()%5)+1.0;
        //    t0=float(randIntCustom(0,4))+1.0;
        //    t0=randf0and1();
            cout<<"rara: ";
        //   t0=abs(float(randf0and1())-0.5);
        //    angle_rotategg = (RandomDoubleX(0.0,1.0)*2.0-1.0)*anggg;
//                angle_rotategg = anggg;
//                t0=RandomDoubleX(0.0,3.0);
        }  */
/*        if(acc_rate>10 && acc_rate<15)
        {
            t0=randf0and1()+2.0;
        } */
    //    if(acc_rate>50) acc_rate=0; 
    /*    else if(acc_rate>50 && acc_rate<60)
        {
//            acc_rate=0;
        //    anggg=float(rand()%90)+120.0;
            anggg=float(randIntCustom(0,60))+120.0;
            float tmp_rv=(2.0*randf0and1()-1.0);
            if(tmp_rv>=0.0)
            {
              angle_rotategg=anggg; // rotate angle  
            }
            else
            {
                angle_rotategg=(-1.0)*anggg; // rotate angle  
            }                
        //    t0=float(rand()%5)+1.0;
            t0=abs(float(randf0and1())-0.5);  
        //    t0=float(randIntCustom(0,4))+1.0;
        //    t0=0.1;              
        }        */    
        float t1=(randf0and1()*2.0-1.0)*t0+0.0;
        float t2=(randf0and1()*2.0-1.0)*t0+0.0;
        float t3=(randf0and1()*2.0-1.0)*t0+0.0;  
    //    cout<<"ang,t0: "<<angle_rotategg<<" ("<<t1<<","<<t2<<","<<t3<<")"<<endl;          
        
        float asin=sin(angle_rotategg);
        float acos=cos(angle_rotategg);
        float u[3][3];
        u[0][0]=acos+awx*awx*(1.0-acos);
        u[0][1]=awx*awy*(1.0-acos)-awz*asin;
        u[0][2]=awx*awz*(1.0-acos)+awy*asin;
        u[1][0]=awx*awy*(1.0-acos)+awz*asin;
        u[1][1]=acos+awy*awy*(1.0-acos);
        u[1][2]=awy*awz*(1.0-acos)-awx*asin;
        u[2][0]=awx*awz*(1.0-acos)-awy*asin;
        u[2][1]=awy*awz*(1.0-acos)+awx*asin;
        u[2][2]=acos+awz*awz*(1.0-acos);            
        // Rotation points
//      int rand_point = rand()%(pnumx);
//      int rand_point = rand()%(pnum);
//      axyz[0]=pointsx[rand_point].ca[0];
//      axyz[1]=pointsx[rand_point].ca[1];
//      axyz[2]=pointsx[rand_point].ca[2];
    //    int rand_point = rand()%6+83 ;
    //    int rand_point = randop ;
    //    axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
    //    axyz[1]=pointsBx[3*rand_point+1].x_[1];
    //    axyz[2]=pointsBx[3*rand_point+1].x_[2];
        axyz[0]= coor_pdb[0];  
        axyz[1]= coor_pdb[1];
        axyz[2]= coor_pdb[2];

    //    beg_p=rand_point; 
//    cout<<"XXX2"<<endl;  

        int tmp3=0;
        vector<poseCoord > fin_mat;
    //    fin_mat = tmp_mat;
    //    if(beg_p==0)
    //    {
        tmp3=0;
        vector<poseCoord>().swap(fin_mat);
        for(int j=0;j<pnum;j++)
        {
            fin_mat.push_back(pointsBx[3*j+0]);
            fin_mat.push_back(pointsBx[3*j+1]);
            fin_mat.push_back(pointsBx[3*j+2]);
            tmp3=tmp3+1;
        }
        old_CC1=new_CC2;   
        old_dE = old_CC1 ;// + old_clash + 0.1*d1/d2;// + olddis_p_m;
    //    old_dE = old_clash + olddis_p_m;                      
        vector<poseCoord > fin_matx;
        fin_matx=fin_mat;
        tmp3=0;
        for(int j=0;j<fin_mat.size();j++) // not last N
        {
            fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
            fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
            fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
            tmp3=tmp3+1;
        }
        vector<poseCoord>().swap(fin_matx);

        for ( int t=0; t<theDensityMapx.density.u1()*theDensityMapx.density.u2()*theDensityMapx.density.u3(); ++t ) {
            rhoC0x[t]=0.0;
            inv_rho_mask0x[t]=1.0;
        }  
        del_ijx=vector<float>(3,0.0);
        atm_jx=vector<float>(3,0.0);
        tmp_i=0; 
    //    string elt_i; 
    //    float bond0max=0.0,bond1max=0.0,bond2max=0.0;
    //    float bond0min=1000.0,bond1min=1000.0,bond2min=1000.0;   
        bond0max=-10000.0;bond1max=-10000.0;bond2max=-10000.0;
        bond0min=10000.0;bond1min=10000.0;bond2min=10000.0;       
        for(int i=0; i<fin_mat.size();i++)
        {
    //        cout<<"YYYY: "<<pointsBx[i].elt_<<endl;;
    //        if(pointsBx[i].elt_ == "CA")
    //        {
            if(fin_mat[i].elt_ == "CA")
            {           
    //            cout<<"RRRR: "<<i;
                vector<float> cartX1;
                vector<float> fracX1;
                elt_i = fin_mat[i].elt_;
                elt_i = elt_i[0];
                OneGaussianScattering sig_j = get_A( elt_i );
                k = sig_j.k( theDensityMapx.effectiveB );
                C = sig_j.C( k );
                if ( C < 1e-6 ) continue;   

                cartX1 = fin_mat[i].x_;           
                MatrixTimesTransVector(theDensityMapx.c2f,cartX1,fracX1);

                // the location of atom in grid ?
                atm_idx[0] = pos_mod (double(fracX1[0]*theDensityMapx.grid[0] - theDensityMapx.origin[0] + 1) , (double)theDensityMapx.grid[0]);
                atm_idx[1] = pos_mod (double(fracX1[1]*theDensityMapx.grid[1] - theDensityMapx.origin[1] + 1) , (double)theDensityMapx.grid[1]);
                atm_idx[2] = pos_mod (double(fracX1[2]*theDensityMapx.grid[2] - theDensityMapx.origin[2] + 1) , (double)theDensityMapx.grid[2]);   
        //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
        //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
        //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
        //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
                for(int z=1;z<=theDensityMapx.density.u3();z++)
                {
                    atm_jx[2] = z;
                    del_ijx[2] =(atm_idx[2]-atm_jx[2])/theDensityMapx.grid[2];
                    if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                    if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                    del_ijx[0] = del_ijx[1] = 0.0;
                    vector<float> frac_tmpz;
                    MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,frac_tmpz);
                    if(square_len(frac_tmpz)> (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;                    
                    for(int y=1;y<=theDensityMapx.density.u2();y++)
                    {
                        atm_jx[1] = y;
                        del_ijx[1] = (atm_idx[1] - atm_jx[1])/theDensityMapx.grid[1];
                        // wrap-around??
                        if ( del_ijx[1] > 0.5 ) del_ijx[1]-=1.0;
                        if ( del_ijx[1] < -0.5 ) del_ijx[1]+=1.0;                        
                        del_ijx[0] = 0.0;
                        vector<float> frac_tmpy;
                        MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,frac_tmpy);
                        if(square_len(frac_tmpy)> (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;                        
                        for(int x=1;x<=theDensityMapx.density.u1();x++)
                        {
                            atm_jx[0] = x;
                            del_ijx[0] = (atm_idx[0] - atm_jx[0])/theDensityMapx.grid[0];
                            // wrap-around??
                            if ( del_ijx[0] > 0.5 ) del_ijx[0]-=1.0;
                            if ( del_ijx[0] < -0.5 ) del_ijx[0]+=1.0;                            
                            vector<float> cart_del_ij2;
                            MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,cart_del_ij2);
                            float d2 = square_len(cart_del_ij2);
                            if(d2 > (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;
                        
                            float atm = C*exp(-k*d2);
    //                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
                            float sigmoid_msk = exp( d2 - (theDensityMapx.ATOM_MASK)*(theDensityMapx.ATOM_MASK)  );
                        //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
                        //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
                            float inv_msk = 1/(1+sigmoid_msk);
                            rhoC0x(x,y,z) += atm;
                            inv_rho_mask0x(x,y,z) *= (1 - inv_msk);

                            if(x>bond0max) bond0max = x;
                            if(y>bond1max) bond1max = y;
                            if(z>bond2max) bond2max = z;
                            if(x<bond0min) bond0min = x;
                            if(y<bond1min) bond1min = y;
                            if(z<bond2min) bond2min = z;                                               

    //                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
    //                                mask2(x,y,z) = 1.0; // problem?
    //                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
    //                                    float atm = C*exp(-k*d2);
    //                                    rhoC(x,y,z) += atm;
    //                                }
    //                            }
                        }
                    }
                } 
    //            cout<<" "<<tmp_i<<" "<<endl;           
                tmp_i = tmp_i + 1;
            }    
    //        }
        } 
        sumC_i2=0.0; sumO_i2=0.0; sumCO_i2=0.0; vol_i2=0.0; CC_i2=0.0;
        sumO2_i2=0.0; sumC2_i2=0.0; varC_i2=0.0; varO_i2=0.0;
        clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
        allpdbmass=0.0;
        allmrcmass=0.0;
        mass_pdb=vector<float>(3,0.0);
        mass_mrc=vector<float>(3,0.0);
        for(int x=int(bond0min);x<=bond0max;x++)
        {
           for(int y=int(bond1min);y<=bond1max;y++)
            {
                for(int z=int(bond2min);z<=bond2max;z++)
                {
                    clc_x2 = rhoC0x(x,y,z);
    //                clc_x2 = tmp_rhc;
                    obs_x2 = theDensityMapx.density(x,y,z);
    //                obs_x2 = tmp_den;
                    eps_x2 = 1-inv_rho_mask0x(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                    sumCO_i2 += eps_x2*clc_x2*obs_x2;
                    sumO_i2  += eps_x2*obs_x2;
                    sumO2_i2 += eps_x2*obs_x2*obs_x2;
                    sumC_i2  += eps_x2*clc_x2;
                    sumC2_i2 += eps_x2*clc_x2*clc_x2;
                    vol_i2   += eps_x2;
                }
            }
        }

        varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
        varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 );
        if ( varC_i2 == 0 || varO_i2 == 0 ) {
            CC_i2 = 0;
        } else {
            CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
        }
        new_CC1=1.0-CC_i2; 
        new_dE = new_CC1 ;                                   
//        cout<<"old_CC1,new_CC1:RRRRRRRRRRRRR "<<old_CC1<<" "<<new_CC1<<endl;
        dE = new_dE-old_dE;
        if(new_dE<old_dE)
        {
            tmp3=0;
            coor_pdb=vector<float>(3,0.0);
            for(int j=0;j<pnum;j++)
            {
                pointsBx[3*j+0]=fin_mat[3*tmp3+0];
                pointsBx[3*j+1]=fin_mat[3*tmp3+1];
                pointsBx[3*j+2]=fin_mat[3*tmp3+2];

                coor_pdb[0]= coor_pdb[0] + fin_mat[3*tmp3+1].x_[0];
                coor_pdb[1]= coor_pdb[1] + fin_mat[3*tmp3+1].x_[1];
                coor_pdb[2]= coor_pdb[2] + fin_mat[3*tmp3+1].x_[2];                    
                tmp3=tmp3+1;
            }
            coor_pdb[0] = coor_pdb[0]/(pnum);
            coor_pdb[1] = coor_pdb[1]/(pnum);
            coor_pdb[2] = coor_pdb[2]/(pnum); 
            if(new_dE<best_E)
            {
                best_E = new_dE;
                vector<poseCoord>().swap(best_model);
                best_model = fin_mat;
            }
        //    cout<<"coor: "<<coor_pdb[0]<<" "<<coor_pdb[1]<<" "<<coor_pdb[2]<<endl;                  
            new_CC2 = new_CC1;
             
//            cout<<"old_CC1,new_CC1: "<<old_dE<<" "<<new_dE<<endl;         
        }else
        {
        //    float tmpx=rand()/double(RAND_MAX);
            float tmpx=randf0and1();
            float mc_v = exp(-dE/(KT));
            if(tmpx< mc_v) // CC must be >0
            {
                tmp3=0;
                coor_pdb=vector<float>(3,0.0);
                for(int j=0;j<pnum;j++)
                {
                    pointsBx[3*j+0]=fin_mat[3*tmp3+0];
                    pointsBx[3*j+1]=fin_mat[3*tmp3+1];
                    pointsBx[3*j+2]=fin_mat[3*tmp3+2];

                    coor_pdb[0]= coor_pdb[0] + fin_mat[3*tmp3+1].x_[0];
                    coor_pdb[1]= coor_pdb[1] + fin_mat[3*tmp3+1].x_[1];
                    coor_pdb[2]= coor_pdb[2] + fin_mat[3*tmp3+1].x_[2];                         
                    tmp3=tmp3+1;
                } 
                coor_pdb[0] = coor_pdb[0]/(pnum);
                coor_pdb[1] = coor_pdb[1]/(pnum);
                coor_pdb[2] = coor_pdb[2]/(pnum); 
            //    cout<<"coor: "<<coor_pdb[0]<<" "<<coor_pdb[1]<<" "<<coor_pdb[2]<<endl;                     
                new_CC2 = new_CC1;
//                cout<<"old_CC1,new_CC1:XX "<<old_dE<<" "<<new_dE<<endl;                    
            }                
        }           
    }
    // rotate 0 and 180
    int tmp3=0;
    coor_pdb=vector<float>(3,0.0);
    for(int j=0;j<pnum;j++)
    {
        coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
        coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
        coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];                         
        tmp3=tmp3+1;
    } 
    coor_pdb[0] = coor_pdb[0]/(pnum);
    coor_pdb[1] = coor_pdb[1]/(pnum);
    coor_pdb[2] = coor_pdb[2]/(pnum); 
    float init_E0=0.0;
    float init_E1=0.0;  
    vector<poseCoord> init_mat0;
    vector<poseCoord> init_mat1;    
    for(int kkk=0;kkk<0;kkk++) // 2
    {
        vector<poseCoord > point_mat;      
        if(kkk==0)
        {
            tmp3=0;
            vector<poseCoord>().swap(point_mat);
            for(int j=0;j<pnum;j++)
            {
                point_mat.push_back(pointsBx[3*j+0]);
                point_mat.push_back(pointsBx[3*j+1]);
                point_mat.push_back(pointsBx[3*j+2]);
                tmp3=tmp3+1;
            }            
        }
        if(kkk==1)
        {
            tmp3=0;
            coor_pdb=vector<float>(3,0.0);
            vector<poseCoord>().swap(point_mat);
            for(int j=0;j<pnum;j++)
            {
                point_mat.push_back(pointsBx[3*j+0]);
                point_mat.push_back(pointsBx[3*j+1]);
                point_mat.push_back(pointsBx[3*j+2]);
                coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
                coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
                coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];                     
                tmp3=tmp3+1;
            }              
            coor_pdb[0] = coor_pdb[0]/(pnum);
            coor_pdb[1] = coor_pdb[1]/(pnum);
            coor_pdb[2] = coor_pdb[2]/(pnum);             
            float asin_theta=2.0*randf0and1()-1.0;
            float acos_theta=sqrt(1.0-asin_theta*asin_theta);
            float apha=2.0*PI*randf0and1();
            float awx=acos_theta*cos(apha);
            float awy=acos_theta*sin(apha);
            float awz=asin_theta;
            // Translation Vector
            float t0=0.0;
            float angle_rotategg=180.0; // rotate angle  
            float t1=(randf0and1()*2.0-1.0)*t0+0.0;
            float t2=(randf0and1()*2.0-1.0)*t0+0.0;
            float t3=(randf0and1()*2.0-1.0)*t0+0.0;  
        //    cout<<"ang,t0: "<<angle_rotategg<<" ("<<t1<<","<<t2<<","<<t3<<")"<<endl;          
            
            float asin=sin(angle_rotategg);
            float acos=cos(angle_rotategg);
            float u[3][3];
            u[0][0]=acos+awx*awx*(1.0-acos);
            u[0][1]=awx*awy*(1.0-acos)-awz*asin;
            u[0][2]=awx*awz*(1.0-acos)+awy*asin;
            u[1][0]=awx*awy*(1.0-acos)+awz*asin;
            u[1][1]=acos+awy*awy*(1.0-acos);
            u[1][2]=awy*awz*(1.0-acos)-awx*asin;
            u[2][0]=awx*awz*(1.0-acos)-awy*asin;
            u[2][1]=awy*awz*(1.0-acos)+awx*asin;
            u[2][2]=acos+awz*awz*(1.0-acos);            
            // Rotation points
    //      int rand_point = rand()%(pnumx);
    //      int rand_point = rand()%(pnum);
    //      axyz[0]=pointsx[rand_point].ca[0];
    //      axyz[1]=pointsx[rand_point].ca[1];
    //      axyz[2]=pointsx[rand_point].ca[2];
        //    int rand_point = rand()%6+83 ;
        //    int rand_point = randop ;
        //    axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
        //    axyz[1]=pointsBx[3*rand_point+1].x_[1];
        //    axyz[2]=pointsBx[3*rand_point+1].x_[2];
            axyz[0]= coor_pdb[0];  
            axyz[1]= coor_pdb[1];
            axyz[2]= coor_pdb[2];

        //    beg_p=rand_point; 
    //    cout<<"XXX2"<<endl;  

        //    vector<poseCoord > fin_mat;
        //    fin_mat = tmp_mat;
        //    if(beg_p==0)
        //    {                   
            vector<poseCoord > fin_matx;
            fin_matx=point_mat;
            tmp3=0;
            for(int j=0;j<point_mat.size();j++) // not last N
            {
                point_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                point_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                point_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                tmp3=tmp3+1;
            }
            vector<poseCoord>().swap(fin_matx);
        } 
        best_E=10.0;
        coor_pdb=vector<float>(3,0.0);   
        for(int j=0;j<pnum;j++)
        {        
            coor_pdb[0]= coor_pdb[0] + point_mat[3*j+1].x_[0];
            coor_pdb[1]= coor_pdb[1] + point_mat[3*j+1].x_[1];
            coor_pdb[2]= coor_pdb[2] + point_mat[3*j+1].x_[2];       
        }
        coor_pdb[0] = coor_pdb[0]/(pnum);
        coor_pdb[1] = coor_pdb[1]/(pnum);
        coor_pdb[2] = coor_pdb[2]/(pnum);      
        for ( int t=0; t<theDensityMapx.density.u1()*theDensityMapx.density.u2()*theDensityMapx.density.u3(); ++t ) {
            rhoC0x[t]=0.0;
            inv_rho_mask0x[t]=1.0;
        }  
        del_ijx=vector<float>(3,0.0);
        atm_jx=vector<float>(3,0.0);
        tmp_i=0; 
    //    string elt_i; 
    //    float bond0max=0.0,bond1max=0.0,bond2max=0.0;
    //    float bond0min=1000.0,bond1min=1000.0,bond2min=1000.0;   
        bond0max=-10000.0;bond1max=-10000.0;bond2max=-10000.0;
        bond0min=10000.0;bond1min=10000.0;bond2min=10000.0;       
        for(int i=0; i<point_mat.size();i++)
        {
    //        cout<<"YYYY: "<<pointsBx[i].elt_<<endl;;
    //        if(pointsBx[i].elt_ == "CA")
    //        {
            if(point_mat[i].elt_ == "CA")
            {           
    //            cout<<"RRRR: "<<i;
                vector<float> cartX1;
                vector<float> fracX1;
                elt_i = point_mat[i].elt_;
                elt_i = elt_i[0];
                OneGaussianScattering sig_j = get_A( elt_i );
                k = sig_j.k( theDensityMapx.effectiveB );
                C = sig_j.C( k );
                if ( C < 1e-6 ) continue;   

                cartX1 = point_mat[i].x_;           
                MatrixTimesTransVector(theDensityMapx.c2f,cartX1,fracX1);

                // the location of atom in grid ?
                atm_idx[0] = pos_mod (double(fracX1[0]*theDensityMapx.grid[0] - theDensityMapx.origin[0] + 1) , (double)theDensityMapx.grid[0]);
                atm_idx[1] = pos_mod (double(fracX1[1]*theDensityMapx.grid[1] - theDensityMapx.origin[1] + 1) , (double)theDensityMapx.grid[1]);
                atm_idx[2] = pos_mod (double(fracX1[2]*theDensityMapx.grid[2] - theDensityMapx.origin[2] + 1) , (double)theDensityMapx.grid[2]);   
        //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
        //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
        //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
        //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
                for(int z=1;z<=theDensityMapx.density.u3();z++)
                {
                    atm_jx[2] = z;
                    del_ijx[2] =(atm_idx[2]-atm_jx[2])/theDensityMapx.grid[2];
                    if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                    if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                    del_ijx[0] = del_ijx[1] = 0.0;
                    vector<float> frac_tmpz;
                    MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,frac_tmpz);
                    if(square_len(frac_tmpz)> (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;                    
                    for(int y=1;y<=theDensityMapx.density.u2();y++)
                    {
                        atm_jx[1] = y;
                        del_ijx[1] = (atm_idx[1] - atm_jx[1])/theDensityMapx.grid[1];
                        // wrap-around??
                        if ( del_ijx[1] > 0.5 ) del_ijx[1]-=1.0;
                        if ( del_ijx[1] < -0.5 ) del_ijx[1]+=1.0;                        
                        del_ijx[0] = 0.0;
                        vector<float> frac_tmpy;
                        MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,frac_tmpy);
                        if(square_len(frac_tmpy)> (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;                        
                        for(int x=1;x<=theDensityMapx.density.u1();x++)
                        {
                            atm_jx[0] = x;
                            del_ijx[0] = (atm_idx[0] - atm_jx[0])/theDensityMapx.grid[0];
                            // wrap-around??
                            if ( del_ijx[0] > 0.5 ) del_ijx[0]-=1.0;
                            if ( del_ijx[0] < -0.5 ) del_ijx[0]+=1.0;                            
                            vector<float> cart_del_ij2;
                            MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,cart_del_ij2);
                            float d2 = square_len(cart_del_ij2);
                            if(d2 > (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;
                        
                            float atm = C*exp(-k*d2);
    //                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
                            float sigmoid_msk = exp( d2 - (theDensityMapx.ATOM_MASK)*(theDensityMapx.ATOM_MASK)  );
                        //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
                        //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
                            float inv_msk = 1/(1+sigmoid_msk);
                            rhoC0x(x,y,z) += atm;
                            inv_rho_mask0x(x,y,z) *= (1 - inv_msk);

                            if(x>bond0max) bond0max = x;
                            if(y>bond1max) bond1max = y;
                            if(z>bond2max) bond2max = z;
                            if(x<bond0min) bond0min = x;
                            if(y<bond1min) bond1min = y;
                            if(z<bond2min) bond2min = z;                                               

    //                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
    //                                mask2(x,y,z) = 1.0; // problem?
    //                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
    //                                    float atm = C*exp(-k*d2);
    //                                    rhoC(x,y,z) += atm;
    //                                }
    //                            }
                        }
                    }
                } 
    //            cout<<" "<<tmp_i<<" "<<endl;           
                tmp_i = tmp_i + 1;
            }    
    //        }
        } 
        sumC_i2=0.0; sumO_i2=0.0; sumCO_i2=0.0; vol_i2=0.0; CC_i2=0.0;
        sumO2_i2=0.0; sumC2_i2=0.0; varC_i2=0.0; varO_i2=0.0;
        clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
    //    allpdbmass=0.0;
    //    allmrcmass=0.0;
    //    mass_pdb=vector<float>(3,0.0);
    //    mass_mrc=vector<float>(3,0.0);
        for(int x=int(bond0min);x<=bond0max;x++)
        {
           for(int y=int(bond1min);y<=bond1max;y++)
            {
                for(int z=int(bond2min);z<=bond2max;z++)
                {
                    clc_x2 = rhoC0x(x,y,z);
    //                clc_x2 = tmp_rhc;
                    obs_x2 = theDensityMapx.density(x,y,z);
    //                obs_x2 = tmp_den;
                    eps_x2 = 1-inv_rho_mask0x(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                    sumCO_i2 += eps_x2*clc_x2*obs_x2;
                    sumO_i2  += eps_x2*obs_x2;
                    sumO2_i2 += eps_x2*obs_x2*obs_x2;
                    sumC_i2  += eps_x2*clc_x2;
                    sumC2_i2 += eps_x2*clc_x2*clc_x2;
                    vol_i2   += eps_x2;
                }
            }
        } 
    /*   for ( int x=0; x<theDensityMapx.density.u1()*theDensityMapx.density.u2()*theDensityMapx.density.u3(); ++x ) {
            clc_x2 = rhoC0x[x];
    //            clc_x2 = tmp_rhc;
            obs_x2 = theDensityMapx.density[x];
    //            obs_x2 = tmp_den;
            eps_x2 = 1-inv_rho_mask0x[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
    //            eps_x2 = 1.0;

            // SMOOTHED
            sumCO_i2 += eps_x2*clc_x2*obs_x2;
            sumO_i2  += eps_x2*obs_x2;
            sumO2_i2 += eps_x2*obs_x2*obs_x2;
            sumC_i2  += eps_x2*clc_x2;
            sumC2_i2 += eps_x2*clc_x2*clc_x2;
            vol_i2   += eps_x2;
        }    */
        varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
        varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 );
    //    cout<<"sumO2_i2: "<<sumO2_i2<<" "<<sumO_i2<<" "<<vol_i2<<endl;
    //    cout<<"varC_i2: "<<varC_i2<<" "<<varO_i2<<endl;
        if ( varC_i2 == 0 || varO_i2 == 0 ) {
            CC_i2 = 0;
        } else {
            CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
        }   
        new_CC2=1.0-CC_i2;
    //    cout<<"new_CC2: "<<new_CC2<<endl;
        for(int iii=0;iii<500;iii++)
        {
    //        cout<<"trdisx,anglex,KT: "<<trdis0n<<" "<<angle0n<<" "<<KT<<endl;
    //        float asin_theta=2.0*RandomDoubleX(0.0,1.0)-1.0;
            float asin_theta=2.0*randf0and1()-1.0;
            float acos_theta=sqrt(1.0-asin_theta*asin_theta);
            float apha=2.0*PI*randf0and1();
            float awx=acos_theta*cos(apha);
            float awy=acos_theta*sin(apha);
            float awz=asin_theta;
            // Translation Vector
            float t0=1.0;
        //    float t0=1.0;

        //    cout<<"t0: "<<t1<<" "<<t2<<" "<<t3<<endl;

            // Rotation matrix
            float anggg=30.0;
            float angle_rotategg=(2.0*randf0and1()-1.0)*anggg; // rotate angle
        //    if(new_CC2>0.95) t0=3.0;
        /*    if(acc_rate>40 && acc_rate<50)
            {
    //            acc_rate=0;
            //    anggg=float(rand()%90)+120.0;
                anggg=float(randIntCustom(0,90))+90.0;
           //     anggg=float(randIntCustom(0,120));
                float tmp_rv=(2.0*randf0and1()-1.0);
                if(tmp_rv>=0.0)
                {
                  angle_rotategg=anggg; // rotate angle  
                }
                else
                {
                    angle_rotategg=(-1.0)*anggg; // rotate angle  
                }
            //    t0=float(rand()%5)+1.0;
            //    t0=float(randIntCustom(0,4))+1.0;
            //    t0=randf0and1();
                cout<<"rara: ";
            //   t0=abs(float(randf0and1())-0.5);
            //    angle_rotategg = (RandomDoubleX(0.0,1.0)*2.0-1.0)*anggg;
    //                angle_rotategg = anggg;
    //                t0=RandomDoubleX(0.0,3.0);
            }  */
    /*        if(acc_rate>10 && acc_rate<15)
            {
                t0=randf0and1()+2.0;
            } */
        //    if(acc_rate>50) acc_rate=0; 
        /*    else if(acc_rate>50 && acc_rate<60)
            {
    //            acc_rate=0;
            //    anggg=float(rand()%90)+120.0;
                anggg=float(randIntCustom(0,60))+120.0;
                float tmp_rv=(2.0*randf0and1()-1.0);
                if(tmp_rv>=0.0)
                {
                  angle_rotategg=anggg; // rotate angle  
                }
                else
                {
                    angle_rotategg=(-1.0)*anggg; // rotate angle  
                }                
            //    t0=float(rand()%5)+1.0;
                t0=abs(float(randf0and1())-0.5);  
            //    t0=float(randIntCustom(0,4))+1.0;
            //    t0=0.1;              
            }        */    
            float t1=(randf0and1()*2.0-1.0)*t0+0.0;
            float t2=(randf0and1()*2.0-1.0)*t0+0.0;
            float t3=(randf0and1()*2.0-1.0)*t0+0.0;  
        //    cout<<"ang,t0: "<<angle_rotategg<<" ("<<t1<<","<<t2<<","<<t3<<")"<<endl;          
            
            float asin=sin(angle_rotategg);
            float acos=cos(angle_rotategg);
            float u[3][3];
            u[0][0]=acos+awx*awx*(1.0-acos);
            u[0][1]=awx*awy*(1.0-acos)-awz*asin;
            u[0][2]=awx*awz*(1.0-acos)+awy*asin;
            u[1][0]=awx*awy*(1.0-acos)+awz*asin;
            u[1][1]=acos+awy*awy*(1.0-acos);
            u[1][2]=awy*awz*(1.0-acos)-awx*asin;
            u[2][0]=awx*awz*(1.0-acos)-awy*asin;
            u[2][1]=awy*awz*(1.0-acos)+awx*asin;
            u[2][2]=acos+awz*awz*(1.0-acos);            
            // Rotation points
    //      int rand_point = rand()%(pnumx);
    //      int rand_point = rand()%(pnum);
    //      axyz[0]=pointsx[rand_point].ca[0];
    //      axyz[1]=pointsx[rand_point].ca[1];
    //      axyz[2]=pointsx[rand_point].ca[2];
        //    int rand_point = rand()%6+83 ;
        //    int rand_point = randop ;
        //    axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
        //    axyz[1]=pointsBx[3*rand_point+1].x_[1];
        //    axyz[2]=pointsBx[3*rand_point+1].x_[2];
            axyz[0]= coor_pdb[0];  
            axyz[1]= coor_pdb[1];
            axyz[2]= coor_pdb[2];

        //    beg_p=rand_point; 
    //    cout<<"XXX2"<<endl;  

            int tmp3=0;
            vector<poseCoord > fin_mat;
        //    fin_mat = tmp_mat;
        //    if(beg_p==0)
        //    {
            tmp3=0;
            vector<poseCoord>().swap(fin_mat);
            for(int j=0;j<pnum;j++)
            {
                fin_mat.push_back(point_mat[3*j+0]);
                fin_mat.push_back(point_mat[3*j+1]);
                fin_mat.push_back(point_mat[3*j+2]);
                tmp3=tmp3+1;
            }
            old_CC1=new_CC2;   
            old_dE = old_CC1 ;// + old_clash + 0.1*d1/d2;// + olddis_p_m;
        //    old_dE = old_clash + olddis_p_m;                      
            vector<poseCoord > fin_matx;
            fin_matx=fin_mat;
            tmp3=0;
            for(int j=0;j<fin_mat.size();j++) // not last N
            {
                fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                tmp3=tmp3+1;
            }
            vector<poseCoord>().swap(fin_matx);

            for ( int t=0; t<theDensityMapx.density.u1()*theDensityMapx.density.u2()*theDensityMapx.density.u3(); ++t ) {
                rhoC0x[t]=0.0;
                inv_rho_mask0x[t]=1.0;
            }  
            del_ijx=vector<float>(3,0.0);
            atm_jx=vector<float>(3,0.0);
            tmp_i=0; 
        //    string elt_i; 
        //    float bond0max=0.0,bond1max=0.0,bond2max=0.0;
        //    float bond0min=1000.0,bond1min=1000.0,bond2min=1000.0;   
            bond0max=-10000.0;bond1max=-10000.0;bond2max=-10000.0;
            bond0min=10000.0;bond1min=10000.0;bond2min=10000.0;       
            for(int i=0; i<fin_mat.size();i++)
            {
        //        cout<<"YYYY: "<<pointsBx[i].elt_<<endl;;
        //        if(pointsBx[i].elt_ == "CA")
        //        {
                if(fin_mat[i].elt_ == "CA")
                {           
        //            cout<<"RRRR: "<<i;
                    vector<float> cartX1;
                    vector<float> fracX1;
                    elt_i = fin_mat[i].elt_;
                    elt_i = elt_i[0];
                    OneGaussianScattering sig_j = get_A( elt_i );
                    k = sig_j.k( theDensityMapx.effectiveB );
                    C = sig_j.C( k );
                    if ( C < 1e-6 ) continue;   

                    cartX1 = fin_mat[i].x_;           
                    MatrixTimesTransVector(theDensityMapx.c2f,cartX1,fracX1);

                    // the location of atom in grid ?
                    atm_idx[0] = pos_mod (double(fracX1[0]*theDensityMapx.grid[0] - theDensityMapx.origin[0] + 1) , (double)theDensityMapx.grid[0]);
                    atm_idx[1] = pos_mod (double(fracX1[1]*theDensityMapx.grid[1] - theDensityMapx.origin[1] + 1) , (double)theDensityMapx.grid[1]);
                    atm_idx[2] = pos_mod (double(fracX1[2]*theDensityMapx.grid[2] - theDensityMapx.origin[2] + 1) , (double)theDensityMapx.grid[2]);   
            //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
            //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
            //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
            //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
                    for(int z=1;z<=theDensityMapx.density.u3();z++)
                    {
                        atm_jx[2] = z;
                        del_ijx[2] =(atm_idx[2]-atm_jx[2])/theDensityMapx.grid[2];
                        if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                        if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                        del_ijx[0] = del_ijx[1] = 0.0;
                        vector<float> frac_tmpz;
                        MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,frac_tmpz);
                        if(square_len(frac_tmpz)> (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;                    
                        for(int y=1;y<=theDensityMapx.density.u2();y++)
                        {
                            atm_jx[1] = y;
                            del_ijx[1] = (atm_idx[1] - atm_jx[1])/theDensityMapx.grid[1];
                            // wrap-around??
                            if ( del_ijx[1] > 0.5 ) del_ijx[1]-=1.0;
                            if ( del_ijx[1] < -0.5 ) del_ijx[1]+=1.0;                        
                            del_ijx[0] = 0.0;
                            vector<float> frac_tmpy;
                            MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,frac_tmpy);
                            if(square_len(frac_tmpy)> (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;                        
                            for(int x=1;x<=theDensityMapx.density.u1();x++)
                            {
                                atm_jx[0] = x;
                                del_ijx[0] = (atm_idx[0] - atm_jx[0])/theDensityMapx.grid[0];
                                // wrap-around??
                                if ( del_ijx[0] > 0.5 ) del_ijx[0]-=1.0;
                                if ( del_ijx[0] < -0.5 ) del_ijx[0]+=1.0;                            
                                vector<float> cart_del_ij2;
                                MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,cart_del_ij2);
                                float d2 = square_len(cart_del_ij2);
                                if(d2 > (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;
                            
                                float atm = C*exp(-k*d2);
        //                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
                                float sigmoid_msk = exp( d2 - (theDensityMapx.ATOM_MASK)*(theDensityMapx.ATOM_MASK)  );
                            //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
                            //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
                                float inv_msk = 1/(1+sigmoid_msk);
                                rhoC0x(x,y,z) += atm;
                                inv_rho_mask0x(x,y,z) *= (1 - inv_msk);

                                if(x>bond0max) bond0max = x;
                                if(y>bond1max) bond1max = y;
                                if(z>bond2max) bond2max = z;
                                if(x<bond0min) bond0min = x;
                                if(y<bond1min) bond1min = y;
                                if(z<bond2min) bond2min = z;                                               

        //                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
        //                                mask2(x,y,z) = 1.0; // problem?
        //                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
        //                                    float atm = C*exp(-k*d2);
        //                                    rhoC(x,y,z) += atm;
        //                                }
        //                            }
                            }
                        }
                    } 
        //            cout<<" "<<tmp_i<<" "<<endl;           
                    tmp_i = tmp_i + 1;
                }    
        //        }
            } 
            sumC_i2=0.0; sumO_i2=0.0; sumCO_i2=0.0; vol_i2=0.0; CC_i2=0.0;
            sumO2_i2=0.0; sumC2_i2=0.0; varC_i2=0.0; varO_i2=0.0;
            clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
            allpdbmass=0.0;
            allmrcmass=0.0;
            mass_pdb=vector<float>(3,0.0);
            mass_mrc=vector<float>(3,0.0);
            for(int x=int(bond0min);x<=bond0max;x++)
            {
               for(int y=int(bond1min);y<=bond1max;y++)
                {
                    for(int z=int(bond2min);z<=bond2max;z++)
                    {
                        clc_x2 = rhoC0x(x,y,z);
        //                clc_x2 = tmp_rhc;
                        obs_x2 = theDensityMapx.density(x,y,z);
        //                obs_x2 = tmp_den;
                        eps_x2 = 1-inv_rho_mask0x(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                        sumCO_i2 += eps_x2*clc_x2*obs_x2;
                        sumO_i2  += eps_x2*obs_x2;
                        sumO2_i2 += eps_x2*obs_x2*obs_x2;
                        sumC_i2  += eps_x2*clc_x2;
                        sumC2_i2 += eps_x2*clc_x2*clc_x2;
                        vol_i2   += eps_x2;
                    }
                }
            }

            varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
            varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 );
            if ( varC_i2 == 0 || varO_i2 == 0 ) {
                CC_i2 = 0;
            } else {
                CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
            }
            new_CC1=1.0-CC_i2; 
            new_dE = new_CC1 ;                                   
//            cout<<"old_CC1,new_CC1:RRRRRRRRRRRRR "<<old_CC1<<" "<<new_CC1<<endl;
            dE = new_dE-old_dE;
            if(new_dE<old_dE)
            {
                tmp3=0;
                coor_pdb=vector<float>(3,0.0);
                for(int j=0;j<pnum;j++)
                {
                    point_mat[3*j+0]=fin_mat[3*tmp3+0];
                    point_mat[3*j+1]=fin_mat[3*tmp3+1];
                    point_mat[3*j+2]=fin_mat[3*tmp3+2];

                    coor_pdb[0]= coor_pdb[0] + fin_mat[3*tmp3+1].x_[0];
                    coor_pdb[1]= coor_pdb[1] + fin_mat[3*tmp3+1].x_[1];
                    coor_pdb[2]= coor_pdb[2] + fin_mat[3*tmp3+1].x_[2];                    
                    tmp3=tmp3+1;
                }
                coor_pdb[0] = coor_pdb[0]/(pnum);
                coor_pdb[1] = coor_pdb[1]/(pnum);
                coor_pdb[2] = coor_pdb[2]/(pnum); 
                if(new_dE<best_E)
                {
                    best_E = new_dE;
                    vector<poseCoord>().swap(best_model);
                    best_model = fin_mat;
                }
            //    cout<<"coor: "<<coor_pdb[0]<<" "<<coor_pdb[1]<<" "<<coor_pdb[2]<<endl;                  
                new_CC2 = new_CC1;
                 
//                cout<<"old_CC1,new_CC1: "<<old_dE<<" "<<new_dE<<endl;         
            }else
            {
            //    float tmpx=rand()/double(RAND_MAX);
                float tmpx=randf0and1();
                float mc_v = exp(-dE/(KT));
                if(tmpx< mc_v) // CC must be >0
                {
                    tmp3=0;
                    coor_pdb=vector<float>(3,0.0);
                    for(int j=0;j<pnum;j++)
                    {
                        point_mat[3*j+0]=fin_mat[3*tmp3+0];
                        point_mat[3*j+1]=fin_mat[3*tmp3+1];
                        point_mat[3*j+2]=fin_mat[3*tmp3+2];

                        coor_pdb[0]= coor_pdb[0] + fin_mat[3*tmp3+1].x_[0];
                        coor_pdb[1]= coor_pdb[1] + fin_mat[3*tmp3+1].x_[1];
                        coor_pdb[2]= coor_pdb[2] + fin_mat[3*tmp3+1].x_[2];                         
                        tmp3=tmp3+1;
                    } 
                    coor_pdb[0] = coor_pdb[0]/(pnum);
                    coor_pdb[1] = coor_pdb[1]/(pnum);
                    coor_pdb[2] = coor_pdb[2]/(pnum); 
                //    cout<<"coor: "<<coor_pdb[0]<<" "<<coor_pdb[1]<<" "<<coor_pdb[2]<<endl;                     
                    new_CC2 = new_CC1;
//                    cout<<"old_CC1,new_CC1:XX "<<old_dE<<" "<<new_dE<<endl;                    
                }                
            }           
        }
        if(kkk==0)
        {
            init_E0=best_E;
            init_mat0 = best_model;
        }
        if(kkk==1)
        {
            init_E1=best_E;
            init_mat1 = best_model;            
        } 
    }
/*    if(init_E0<init_E1)
    {
        vector<poseCoord>().swap(pointsBx);
        pointsBx= init_mat0;
    } else
    {
        vector<poseCoord>().swap(pointsBx);
        pointsBx= init_mat1;
    }  */
//    init0max=bond0max;init1max=bond1max;init2max=bond2max;
//    init0min=bond0min;init1min=bond1min;init2min=bond2min;     
/*    for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
//        cout<<theDensityMap.density[x]<<" ";
        if(theDensityMap.density[x]<(3.0/5.0)*max_den) theDensityMap.density[x]=0; 
    }    */
/*    for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
//        cout<<theDensityMap.density[x]<<" ";
        if(rhoC0[x]<(3.0/5.0)*max_rhc) rhoC0[x]=0; 
    }     */

   
    // ss information
/*    vector<string> ss_inf0(pnum);
    vector<vector<float> > CA_cor0(pnum,vector<float>(3,0.0));
    cout<<"RRRRRzz"<<endl;
    tmp_i=0;
    for(int i=0; i<pointsBx.size();i++)
    {
//        cout<<CA_cor0[tmp_i][0]<<'\t'<<CA_cor0[tmp_i][1]<<'\t'<<CA_cor0[tmp_i][2]<<endl;
        if(pointsBx[i].elt_ == "CA")
        {
            CA_cor0[tmp_i]=pointsBx[i].x_;
//            cout<<CA_cor0[tmp_i][0]<<'\t'<<CA_cor0[tmp_i][1]<<'\t'<<CA_cor0[tmp_i][2]<<endl;
            tmp_i=tmp_i+1;
//            cout<<pointsBx[i].x_[0]<<'\t'<<pointsBx[i].x_[1]<<'\t'<<pointsBx[i].x_[2]<<endl;
        }
    }
    make_secx(CA_cor0,pnum,ss_inf0);  */
    vector<boneinfo> bb(pnum);
    vector<boneinfo> bbh0(pnum);
    vector<poseCoord> proseq;
    proseq=pointsBx;
    int numbb=pnum;
    int tmp_tt=0;
    for(int j=0;j<pnum;j++)
    {
        bb[tmp_tt].indn = 3*tmp_tt+0;
        bb[tmp_tt].indca = 3*tmp_tt+1;
        bb[tmp_tt].indc = 3*tmp_tt+2; 
        string resng = pointsrenm[tmp_tt];
        char resny = resng[0];
        bb[tmp_tt].resid = aminoid(resny);        
        tmp_tt=tmp_tt+1;
    }     
    cout<<"xxExxx"<<endl;
    calcssennhoc(bb,pnum,pointsBx); 

/*    for(int i=0;i<pnum;i++)
    {
//        cout<<numtoss(bb[i].sst)<<" ";
        cout<<":"<<pointsBx[3*i+0].elt_<<" "<<pointsBx[3*i+1].elt_<<" "<<pointsBx[3*i+2].elt_<<endl;
    }    */
    cout<<endl;
/*    for(int i=0;i<pnum;i++)
    {
        if(i==0) continue;
        int tyx=i-1;
        int tyy=i;
        int tyz=i+1;
        if(tyz>=pnum-1) break;
        if(bb[tyy].sst!=bb.sst[tyx] && bb[tyy].sst!=bb[tyz].sst)  bb[i].sst='C';
    } */
    for(int i=0;i<pnum;i++)
    {
        cout<<i<<":"<<numtoss(bb[i].sst)<<" ";
    }    
    bbh0=bb;
    cout<<endl;    

    // domain align
    // rotate around point
//    int beg_p=dm_beg;
//    int end_p=dm_end;
    int randop=-1; // rotate point
    int randopx=-1;
//    if(beg_p>3) randop=beg_p;
//    if(end_p<(pnum-4) && end_p>4) randopx=end_p;
    cout<<"CCCC"<<endl; 

    ObjexxFCL::FArray3D< float > rhoC01;
    ObjexxFCL::FArray3D< float > inv_rho_mask1;
    rhoC01.dimension(theDensityMapx.density.u1() , theDensityMapx.density.u2() , theDensityMapx.density.u3());
    inv_rho_mask1.dimension(theDensityMapx.density.u1() , theDensityMapx.density.u2() , theDensityMapx.density.u3());
    for ( int t=0; t<theDensityMapx.density.u1()*theDensityMapx.density.u2()*theDensityMapx.density.u3(); ++t ) {
        rhoC01[t]=0.0;
        inv_rho_mask1[t]=1.0;
    }

    allpdbmass=0.0;
    allmrcmass=0.0;
    float olddis_p_m=0.0;
    float newdis_p_m=0.0;
    int acc_rate=0; // not acc_rate;
/*    int large_dom0= 0;
    int large_dom1= 0;
    int start_dom = 0;
    int med_dom = 0;
    int end_dom = 0;
    if(dm_beg<3)
    {
        start_dom = 0;
        med_dom = dm_end;
        end_dom = pnum;
    } else
    {
        start_dom = 0;
        med_dom = dm_beg;
        end_dom = dm_end;
    }
    int first_run_dom_start = 0;
    int first_run_dom_end = 0;
    int second_run_dom_start = 0;
    int second_run_dom_end = 0;        
    if((med_dom-start_dom)>=(end_dom-med_dom))
    {
        first_run_dom_start = start_dom;
        first_run_dom_end = med_dom;
        second_run_dom_start = med_dom;
        second_run_dom_end = end_dom;
    } else
    {
        first_run_dom_start = med_dom;
        first_run_dom_end = end_dom;
        second_run_dom_start = start_dom;
        second_run_dom_end = med_dom;            
    }    */
    for(int tttt=0;tttt<0;tttt++)// 2
    {

    //    beg_p=0;
    //    end_p=pnum;
        ObjexxFCL::FArray3D< float > rhoC01x;
        ObjexxFCL::FArray3D< float > inv_rho_mask1x;
        rhoC01x.dimension(theDensityMap.density.u1() , theDensityMap.density.u2() , theDensityMap.density.u3());
        inv_rho_mask1x.dimension(theDensityMap.density.u1() , theDensityMap.density.u2() , theDensityMap.density.u3());
        for ( int t=0; t<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++t ) {
            rhoC01x[t]=0.0;
            inv_rho_mask1x[t]=1.0;
        }
        vector<poseCoord > dom_mat;
        coor_pdb=vector<float>(3,0.0);
        for(int j=0;j<pnum;j++)
        {
            dom_mat.push_back(pointsBx[3*j+0]);
            dom_mat.push_back(pointsBx[3*j+1]);
            dom_mat.push_back(pointsBx[3*j+2]);
            coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
            coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
            coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];
        }
        coor_pdb[0] = coor_pdb[0]/(pnum);
        coor_pdb[1] = coor_pdb[1]/(pnum);
        coor_pdb[2] = coor_pdb[2]/(pnum);

        del_ijx=vector<float>(3,0.0);
        atm_jx=vector<float>(3,0.0);           
    //    vector<vector<float> >().swap(atm_idx);
        atm_idx=vector<float>(3,0.0);
        bond0max=-10000.0;bond1max=-10000.0;bond2max=-10000.0;
        bond0min=10000.0;bond1min=10000.0;bond2min=10000.0;     
        tmp_i=0;       
        for ( int t=0; t<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++t ) {
            rhoC01x[t]=0.0;
            inv_rho_mask1x[t]=1.0;
        }        
        for(int i=0;i<dom_mat.size();i++)  // rotate fragment from beg_p to end_p
        {
             if(dom_mat[i].elt_ == "CA")
            {           
    //            cout<<"RRRR: "<<i;
                vector<float> cartX1;
                vector<float> fracX1;
                elt_i = dom_mat[i].elt_;
                elt_i = elt_i[0];
                OneGaussianScattering sig_j = get_A( elt_i );
                k = sig_j.k( theDensityMap.effectiveB );
                C = sig_j.C( k );
                if ( C < 1e-6 ) continue;   

                cartX1 = dom_mat[i].x_;           
                MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);

                // the location of atom in grid ?
                atm_idx[0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                atm_idx[1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                atm_idx[2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);   
        //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
        //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
        //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
        //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
                for(int z=1;z<=theDensityMap.density.u3();z++)
                {
                    atm_jx[2] = z;
                    del_ijx[2] =(atm_idx[2]-atm_jx[2])/theDensityMap.grid[2];
                    if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                    if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                    del_ijx[0] = del_ijx[1] = 0.0;
                    vector<float> frac_tmpz;
                    MatrixTimesTransVector(theDensityMap.f2c,del_ijx,frac_tmpz);
                    if(square_len(frac_tmpz)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                    
                    for(int y=1;y<=theDensityMap.density.u2();y++)
                    {
                        atm_jx[1] = y;
                        del_ijx[1] = (atm_idx[1] - atm_jx[1])/theDensityMap.grid[1];
                        // wrap-around??
                        if ( del_ijx[1] > 0.5 ) del_ijx[1]-=1.0;
                        if ( del_ijx[1] < -0.5 ) del_ijx[1]+=1.0;                        
                        del_ijx[0] = 0.0;
                        vector<float> frac_tmpy;
                        MatrixTimesTransVector(theDensityMap.f2c,del_ijx,frac_tmpy);
                        if(square_len(frac_tmpy)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                        
                        for(int x=1;x<=theDensityMap.density.u1();x++)
                        {
                            atm_jx[0] = x;
                            del_ijx[0] = (atm_idx[0] - atm_jx[0])/theDensityMap.grid[0];
                            // wrap-around??
                            if ( del_ijx[0] > 0.5 ) del_ijx[0]-=1.0;
                            if ( del_ijx[0] < -0.5 ) del_ijx[0]+=1.0;                            
                            vector<float> cart_del_ij2;
                            MatrixTimesTransVector(theDensityMap.f2c,del_ijx,cart_del_ij2);
                            float d2 = square_len(cart_del_ij2);
                            if(d2 > (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;
                        
                            float atm = C*exp(-k*d2);
    //                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
                            float sigmoid_msk = exp( d2 - (theDensityMap.ATOM_MASK)*(theDensityMap.ATOM_MASK)  );
                        //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
                        //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
                            float inv_msk = 1/(1+sigmoid_msk);
                            rhoC01x(x,y,z) += atm;
                            inv_rho_mask1x(x,y,z) *= (1 - inv_msk);


                            if(x>bond0max) bond0max = x;
                            if(y>bond1max) bond1max = y;
                            if(z>bond2max) bond2max = z;
                            if(x<bond0min) bond0min = x;
                            if(y<bond1min) bond1min = y;
                            if(z<bond2min) bond2min = z;
    //                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
    //                                mask2(x,y,z) = 1.0; // problem?
    //                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
    //                                    float atm = C*exp(-k*d2);
    //                                    rhoC(x,y,z) += atm;
    //                                }
    //                            }
                        }
                    }
                } 
    //            cout<<" "<<tmp_i<<" "<<endl;           
                tmp_i = tmp_i + 1;
            } 
        } 
        sumC_i2=0.0; sumO_i2=0.0; sumCO_i2=0.0; vol_i2=0.0; CC_i2=0.0;
        sumO2_i2=0.0; sumC2_i2=0.0; varC_i2=0.0; varO_i2=0.0;
        clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
        allpdbmass=0.0;
        allmrcmass=0.0;
        mass_pdb=vector<float>(3,0.0);
        mass_mrc=vector<float>(3,0.0);
        for(int x=int(bond0min);x<=bond0max;x++)
        {
           for(int y=int(bond1min);y<=bond1max;y++)
            {
                for(int z=int(bond2min);z<=bond2max;z++)
                {
                    clc_x2 = rhoC01x(x,y,z);
    //                clc_x2 = tmp_rhc;
                    obs_x2 = theDensityMap.density(x,y,z);
    //                obs_x2 = tmp_den;
                    eps_x2 = 1-inv_rho_mask1x(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                    sumCO_i2 += eps_x2*clc_x2*obs_x2;
                    sumO_i2  += eps_x2*obs_x2;
                    sumO2_i2 += eps_x2*obs_x2*obs_x2;
                    sumC_i2  += eps_x2*clc_x2;
                    sumC2_i2 += eps_x2*clc_x2*clc_x2;
                    vol_i2   += eps_x2;
                    mass_pdb[0]=  mass_pdb[0] + eps_x2*clc_x2*x;
                    mass_pdb[1]=  mass_pdb[1] + eps_x2*clc_x2*y;
                    mass_pdb[2]=  mass_pdb[2] + eps_x2*clc_x2*z;
//                    mass_mrc[0]=  mass_mrc[0] + eps_x2*obs_x2*x;
//                    mass_mrc[1]=  mass_mrc[1] + eps_x2*obs_x2*y;
//                    mass_mrc[2]=  mass_mrc[2] + eps_x2*obs_x2*z;
                    allpdbmass = allpdbmass + eps_x2*clc_x2;
//                    allmrcmass = allmrcmass + eps_x2*obs_x2;
                }
            }
        }
        mass_pdb[0] = mass_pdb[0]/allpdbmass;
        mass_pdb[1] = mass_pdb[1]/allpdbmass;
        mass_pdb[2] = mass_pdb[2]/allpdbmass;
//        mass_mrc[0] = mass_mrc[0]/allmrcmass;
//        mass_mrc[1] = mass_mrc[1]/allmrcmass;
//        mass_mrc[2] = mass_mrc[2]/allmrcmass;
    //    olddis_p_m=Distance_point(mass_pdb,mass_mrc);   
        /*            for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
            //            tmp_den=0.0;
            //            tmp_rhc =0.0;
            //            if(theDensityMap.density[x] > (2.0/3.0)*max_den) tmp_den= theDensityMap.density[x];
            //            if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
            //            if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                        // fetch this point
                        clc_x2 = rhoC2[x];
            //            clc_x2 = tmp_rhc;
                        obs_x2 = theDensityMap.density[x];
            //            obs_x2 = tmp_den;
                        eps_x2 = 1-inv_rho_mask2[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
            //            eps_x2 = 1.0;

                        // SMOOTHED
                        sumCO_i2 += eps_x2*clc_x2*obs_x2;
                        sumO_i2  += eps_x2*obs_x2;
                        sumO2_i2 += eps_x2*obs_x2*obs_x2;
                        sumC_i2  += eps_x2*clc_x2;
                        sumC2_i2 += eps_x2*clc_x2*clc_x2;
                        vol_i2   += eps_x2;
                    }
                    varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
                    varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
                    if ( varC_i2 == 0 || varO_i2 == 0 ) {
                        CC_i2 = 0;
                    } else {
                        CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
                    }   */           
        varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
        varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 );
        if ( varC_i2 == 0 || varO_i2 == 0 ) {
            CC_i2 = 0;
        } else {
            CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
        }   
        new_CC2=1.0-CC_i2;

//        float Tend[5]={0.01,0.005,0.002,0.008,0.001};
//        float Pend[5]={1.0,3.0,0.9,1.0,0.5};
//        float Aend[5]={30.0,15.0,120.0,30.0,10.0};     
//        float Tend[5]={0.01,0.005,0.002,0.002,0.001};
//        float Pend[5]={1.0,0.9,1.0,3.0,0.5};
//        float Aend[5]={30.0,10.0,120.0,120.0,10.0};
//        int Rend[5]={600,300,600,600,300};     
//        float Tend[7]={0.01,0.005,0.002,0.002,0.001,0.05,0.001};
//        float Pend[7]={1.0,0.9,1.0,3.0,0.8,1.0,0.8};
//        float Aend[7]={30.0,10.0,120.0,120.0,10.0,60.0,10.0};
//        int Rend[7]={600,300,600,600,300,300,600};   
//        float Tend[5]={0.01,0.005,0.008,0.002,0.001};
//        float Pend[5]={3.0,1.0,1.65,0.6,0.5};
//        float Aend[5]={30.0,10.0,30.0,90.0,10.0};
//        float Tend[5]={0.01,0.002,0.008,0.005,0.001};
//        float Pend[5]={1.0,3.0,0.9,1.6,0.5};
//        float Aend[5]={30.0,10.0,15.0,90.0,10.0};        
    vector<poseCoord>().swap(best_model);
    best_E=10.0;        
    for(int ttt=0;ttt<7;ttt++) // 7
    {
        float KT0x =0.1;
        float KT0y =0.01;
        float KT00= pow(float(KT0y/KT0x),float(float(ttt)/float(6)));
        float KT0=KT0x*float(KT00);

        float KT1x =0.01;
        float KT1y =0.001;
        float KT10= pow(float(KT1y/KT1x),float(float(ttt)/float(6)));
        float KT1=KT1x*float(KT10);


    //    KT=Tend[ttt];
    //    cout<<"KT: "<<KT<<endl;
    //    KT=0.005;

//        float angle0 =180.0;
    //    float angle0 =30.0;
//        float anglen =10;
    //    float anglen =5.0;
    //    float anglex= pow(float(anglen/angle0),float(float(ttt)/float(4)));
    //    float angle0n=angle0*float(anglex);
    //    float angle0n=Aend[ttt];
    //    cout<<"A: "<<angle0n<<endl;
    //    float angle0n=30.0;

   //     float trdis0 =3.0;
   //     float trdisn=0.5;
   //     float trdisx= pow(float(trdisn/trdis0),float(float(ttt)/float(4)));
    //    float trdis0n=trdis0*float(trdisx);
    //    float trdis0n=Pend[ttt];
    //    cout<<"T: "<<trdis0n<<endl;
    //    float trdis0n=2.0;      

    //    cout<<"trdisx,anglex,KT 0: "<<trdis0n<<" "<<angle0n<<" "<<KT<<endl;  

        float dis_p=0.0;
        bool tp_clashx= false;
        old_vwd =0.0;
        old_clash = 0.0;
        old_dE =0.0;
        new_dE =0.0;
        axyz=vector<float>(3,0.0);  
        dE=0.0;  
        bond0max=-10000.0;bond1max=-10000.0;bond2max=-10000.0;
        bond0min=10000.0;bond1min=10000.0;bond2min=10000.0;            
    //    vector<vector<float> > fin_matx;
        new_CC1=0.0;
        old_CC1=0.0;
        new_CC2=0.0;
        new_dE2=0.0;
        tmp_i=0;     
        for ( int t=0; t<theDensityMapx.density.u1()*theDensityMapx.density.u2()*theDensityMapx.density.u3(); ++t ) {
            rhoC01[t]=0.0;
            inv_rho_mask1[t]=1.0;
        }   

        vector<poseCoord > fin_maty;
        coor_pdb=vector<float>(3,0.0);
        for(int j=0;j<pnum;j++)
        {
            fin_maty.push_back(pointsBx[3*j+0]);
            fin_maty.push_back(pointsBx[3*j+1]);
            fin_maty.push_back(pointsBx[3*j+2]);
            coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
            coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
            coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];
        }
        best_model = fin_maty;
        coor_pdb[0] = coor_pdb[0]/(pnum);
        coor_pdb[1] = coor_pdb[1]/(pnum);
        coor_pdb[2] = coor_pdb[2]/(pnum);
        cout<<"coor: "<<coor_pdb[0]<<" "<<coor_pdb[1]<<" "<<coor_pdb[2]<<endl;
        del_ijx=vector<float>(3,0.0);
        atm_jx=vector<float>(3,0.0);         
//        vector<vector<float> >().swap(atm_idx);
        atm_idx=vector<float>(3,0.0);
        for(int i=0;i<fin_maty.size();i++)  // rotate fragment from beg_p to end_p
        {
             if(fin_maty[i].elt_ == "CA")
            {           
    //            cout<<"RRRR: "<<i;
                vector<float> cartX1;
                vector<float> fracX1;
                elt_i = fin_maty[i].elt_;
                elt_i = elt_i[0];
                OneGaussianScattering sig_j = get_A( elt_i );
                k = sig_j.k( theDensityMapx.effectiveB );
                C = sig_j.C( k );
                if ( C < 1e-6 ) continue;   

                cartX1 = fin_maty[i].x_;           
                MatrixTimesTransVector(theDensityMapx.c2f,cartX1,fracX1);

                // the location of atom in grid ?
                atm_idx[0] = pos_mod (double(fracX1[0]*theDensityMapx.grid[0] - theDensityMapx.origin[0] + 1) , (double)theDensityMapx.grid[0]);
                atm_idx[1] = pos_mod (double(fracX1[1]*theDensityMapx.grid[1] - theDensityMapx.origin[1] + 1) , (double)theDensityMapx.grid[1]);
                atm_idx[2] = pos_mod (double(fracX1[2]*theDensityMapx.grid[2] - theDensityMapx.origin[2] + 1) , (double)theDensityMapx.grid[2]);   
        //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
        //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
        //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
        //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
                for(int z=1;z<=theDensityMapx.density.u3();z++)
                {
                    atm_jx[2] = z;
                    del_ijx[2] =(atm_idx[2]-atm_jx[2])/theDensityMapx.grid[2];
                    if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                    if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                    del_ijx[0] = del_ijx[1] = 0.0;
                    vector<float> frac_tmpz;
                    MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,frac_tmpz);
                    if(square_len(frac_tmpz)> (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;                    
                    for(int y=1;y<=theDensityMapx.density.u2();y++)
                    {
                        atm_jx[1] = y;
                        del_ijx[1] = (atm_idx[1] - atm_jx[1])/theDensityMapx.grid[1];
                        // wrap-around??
                        if ( del_ijx[1] > 0.5 ) del_ijx[1]-=1.0;
                        if ( del_ijx[1] < -0.5 ) del_ijx[1]+=1.0;                        
                        del_ijx[0] = 0.0;
                        vector<float> frac_tmpy;
                        MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,frac_tmpy);
                        if(square_len(frac_tmpy)> (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;                        
                        for(int x=1;x<=theDensityMapx.density.u1();x++)
                        {
                            atm_jx[0] = x;
                            del_ijx[0] = (atm_idx[0] - atm_jx[0])/theDensityMapx.grid[0];
                            // wrap-around??
                            if ( del_ijx[0] > 0.5 ) del_ijx[0]-=1.0;
                            if ( del_ijx[0] < -0.5 ) del_ijx[0]+=1.0;                            
                            vector<float> cart_del_ij2;
                            MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,cart_del_ij2);
                            float d2 = square_len(cart_del_ij2);
                            if(d2 > (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;
                        
                            float atm = C*exp(-k*d2);
    //                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
                            float sigmoid_msk = exp( d2 - (theDensityMapx.ATOM_MASK)*(theDensityMapx.ATOM_MASK)  );
                        //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
                        //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
                            float inv_msk = 1/(1+sigmoid_msk);
                            rhoC01(x,y,z) += atm;
                            inv_rho_mask1(x,y,z) *= (1 - inv_msk);


                            if(x>bond0max) bond0max = x;
                            if(y>bond1max) bond1max = y;
                            if(z>bond2max) bond2max = z;
                            if(x<bond0min) bond0min = x;
                            if(y<bond1min) bond1min = y;
                            if(z<bond2min) bond2min = z;
    //                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
    //                                mask2(x,y,z) = 1.0; // problem?
    //                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
    //                                    float atm = C*exp(-k*d2);
    //                                    rhoC(x,y,z) += atm;
    //                                }
    //                            }
                        }
                    }
                } 
    //            cout<<" "<<tmp_i<<" "<<endl;           
                tmp_i = tmp_i + 1;
            } 
        } 

        // bandlimit mask at 'radius'
    /*    ObjexxFCL::FArray3D< std::complex<double> > Finv_rho_mask1;
        fourier::fft3(inv_rho_mask1, Finv_rho_mask1);
        //int H,K,L;
        for ( int z=1; z<=(int)theDensityMap.grid[2]; ++z ) {
            int H = (z < (int)theDensityMap.grid[2]/2) ? z-1 : z-theDensityMap.grid[2] - 1;
            for ( int y=1; y<=(int)theDensityMap.grid[1]; ++y ) {
                int K = (y < (int)theDensityMap.grid[1]/2) ? y-1 : y-theDensityMap.grid[1]-1;
                for ( int x=1; x<=(int)theDensityMap.grid[0]; ++x ) {
                    int L = (x < (int)theDensityMap.grid[0]/2) ? x-1 : x-theDensityMap.grid[0]-1;
                    float S2c =  theDensityMap.S2(H,K,L);

                    // exp fade
                    float highres_limit=500.0;
                    float scale = exp(-S2c*(highres_limit*highres_limit));
                    Finv_rho_mask1(x,y,z) *= scale;
                }
            }
        }
        fourier::ifft3(Finv_rho_mask1, inv_rho_mask1);      */
          
        sumC_i2=0.0; sumO_i2=0.0; sumCO_i2=0.0; vol_i2=0.0; CC_i2=0.0;
        sumO2_i2=0.0; sumC2_i2=0.0; varC_i2=0.0; varO_i2=0.0;
        clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
        allpdbmass=0.0;
        allmrcmass=0.0;
        mass_pdb=vector<float>(3,0.0);
        mass_mrc=vector<float>(3,0.0);
        for(int x=int(bond0min);x<=bond0max;x++)
        {
           for(int y=int(bond1min);y<=bond1max;y++)
            {
                for(int z=int(bond2min);z<=bond2max;z++)
                {
                    clc_x2 = rhoC01(x,y,z);
    //                clc_x2 = tmp_rhc;
                    obs_x2 = theDensityMapx.density(x,y,z);
    //                obs_x2 = tmp_den;
                    eps_x2 = 1-inv_rho_mask1(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                    sumCO_i2 += eps_x2*clc_x2*obs_x2;
                    sumO_i2  += eps_x2*obs_x2;
                    sumO2_i2 += eps_x2*obs_x2*obs_x2;
                    sumC_i2  += eps_x2*clc_x2;
                    sumC2_i2 += eps_x2*clc_x2*clc_x2;
                    vol_i2   += eps_x2;
                    mass_pdb[0]=  mass_pdb[0] + eps_x2*clc_x2*x;
                    mass_pdb[1]=  mass_pdb[1] + eps_x2*clc_x2*y;
                    mass_pdb[2]=  mass_pdb[2] + eps_x2*clc_x2*z;
                    mass_mrc[0]=  mass_mrc[0] + eps_x2*obs_x2*x;
                    mass_mrc[1]=  mass_mrc[1] + eps_x2*obs_x2*y;
                    mass_mrc[2]=  mass_mrc[2] + eps_x2*obs_x2*z;
                    allpdbmass = allpdbmass + eps_x2*clc_x2;
                    allmrcmass = allmrcmass + eps_x2*obs_x2;
                }
            }
        }
        mass_pdb[0] = mass_pdb[0]/allpdbmass;
        mass_pdb[1] = mass_pdb[1]/allpdbmass;
        mass_pdb[2] = mass_pdb[2]/allpdbmass;
        mass_mrc[0] = mass_mrc[0]/allmrcmass;
        mass_mrc[1] = mass_mrc[1]/allmrcmass;
        mass_mrc[2] = mass_mrc[2]/allmrcmass;
    //    olddis_p_m=Distance_point(mass_pdb,mass_mrc);   
        /*            for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
            //            tmp_den=0.0;
            //            tmp_rhc =0.0;
            //            if(theDensityMap.density[x] > (2.0/3.0)*max_den) tmp_den= theDensityMap.density[x];
            //            if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
            //            if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                        // fetch this point
                        clc_x2 = rhoC2[x];
            //            clc_x2 = tmp_rhc;
                        obs_x2 = theDensityMap.density[x];
            //            obs_x2 = tmp_den;
                        eps_x2 = 1-inv_rho_mask2[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
            //            eps_x2 = 1.0;

                        // SMOOTHED
                        sumCO_i2 += eps_x2*clc_x2*obs_x2;
                        sumO_i2  += eps_x2*obs_x2;
                        sumO2_i2 += eps_x2*obs_x2*obs_x2;
                        sumC_i2  += eps_x2*clc_x2;
                        sumC2_i2 += eps_x2*clc_x2*clc_x2;
                        vol_i2   += eps_x2;
                    }
                    varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
                    varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
                    if ( varC_i2 == 0 || varO_i2 == 0 ) {
                        CC_i2 = 0;
                    } else {
                        CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
                    }   */           
        varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
        varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 );
        if ( varC_i2 == 0 || varO_i2 == 0 ) {
            CC_i2 = 0;
        } else {
            CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
        }   
        new_CC2=1.0-CC_i2;   
    //    acc_rate=0;
    //    cout<<"Rend: "<<Rend[ttt]<<endl;
        for(int iii=0;iii<500;iii++)
        {
    //        cout<<"trdisx,anglex,KT: "<<trdis0n<<" "<<angle0n<<" "<<KT<<endl;
    //        float asin_theta=2.0*RandomDoubleX(0.0,1.0)-1.0;
//            float KT0 =0.1;
//            float KTn =0.001;
            float KTx= pow(float(KT1/KT0),float(float(iii)/float(499)));
            KT=KT0*float(KTx); 

            float asin_theta=2.0*randf0and1()-1.0;
            float acos_theta=sqrt(1.0-asin_theta*asin_theta);
            float apha=2.0*PI*randf0and1();
            float awx=acos_theta*cos(apha);
            float awy=acos_theta*sin(apha);
            float awz=asin_theta;
            // Translation Vector
            float t0=1.0;
        //    float t0=1.0;

        //    cout<<"t0: "<<t1<<" "<<t2<<" "<<t3<<endl;

            // Rotation matrix
            float anggg=30.0;
            float angle_rotategg=(2.0*randf0and1()-1.0)*anggg; // rotate angle
        //    if(new_CC2>0.95) t0=3.0;
        /*    if(acc_rate>40 && acc_rate<50)
            {
    //            acc_rate=0;
            //    anggg=float(rand()%90)+120.0;
                anggg=float(randIntCustom(0,90))+90.0;
           //     anggg=float(randIntCustom(0,120));
                float tmp_rv=(2.0*randf0and1()-1.0);
                if(tmp_rv>=0.0)
                {
                  angle_rotategg=anggg; // rotate angle  
                }
                else
                {
                    angle_rotategg=(-1.0)*anggg; // rotate angle  
                }
            //    t0=float(rand()%5)+1.0;
            //    t0=float(randIntCustom(0,4))+1.0;
            //    t0=randf0and1();
                cout<<"rara: ";
            //   t0=abs(float(randf0and1())-0.5);
            //    angle_rotategg = (RandomDoubleX(0.0,1.0)*2.0-1.0)*anggg;
//                angle_rotategg = anggg;
//                t0=RandomDoubleX(0.0,3.0);
            }  */
    /*        if(acc_rate>10 && acc_rate<15)
            {
                t0=randf0and1()+2.0;
            } */
        //    if(acc_rate>50) acc_rate=0; 
        /*    else if(acc_rate>50 && acc_rate<60)
            {
    //            acc_rate=0;
            //    anggg=float(rand()%90)+120.0;
                anggg=float(randIntCustom(0,60))+120.0;
                float tmp_rv=(2.0*randf0and1()-1.0);
                if(tmp_rv>=0.0)
                {
                  angle_rotategg=anggg; // rotate angle  
                }
                else
                {
                    angle_rotategg=(-1.0)*anggg; // rotate angle  
                }                
            //    t0=float(rand()%5)+1.0;
                t0=abs(float(randf0and1())-0.5);  
            //    t0=float(randIntCustom(0,4))+1.0;
            //    t0=0.1;              
            }        */    
            float t1=(randf0and1()*2.0-1.0)*t0+0.0;
            float t2=(randf0and1()*2.0-1.0)*t0+0.0;
            float t3=(randf0and1()*2.0-1.0)*t0+0.0;  
        //    cout<<"ang,t0: "<<angle_rotategg<<" ("<<t1<<","<<t2<<","<<t3<<")"<<endl;          
            
            float asin=sin(angle_rotategg);
            float acos=cos(angle_rotategg);
            float u[3][3];
            u[0][0]=acos+awx*awx*(1.0-acos);
            u[0][1]=awx*awy*(1.0-acos)-awz*asin;
            u[0][2]=awx*awz*(1.0-acos)+awy*asin;
            u[1][0]=awx*awy*(1.0-acos)+awz*asin;
            u[1][1]=acos+awy*awy*(1.0-acos);
            u[1][2]=awy*awz*(1.0-acos)-awx*asin;
            u[2][0]=awx*awz*(1.0-acos)-awy*asin;
            u[2][1]=awy*awz*(1.0-acos)+awx*asin;
            u[2][2]=acos+awz*awz*(1.0-acos);            
            // Rotation points
    //      int rand_point = rand()%(pnumx);
    //      int rand_point = rand()%(pnum);
    //      axyz[0]=pointsx[rand_point].ca[0];
    //      axyz[1]=pointsx[rand_point].ca[1];
    //      axyz[2]=pointsx[rand_point].ca[2];
        //    int rand_point = rand()%6+83 ;
        //    int rand_point = randop ;
        //    axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
        //    axyz[1]=pointsBx[3*rand_point+1].x_[1];
        //    axyz[2]=pointsBx[3*rand_point+1].x_[2];
            axyz[0]= coor_pdb[0];  
            axyz[1]= coor_pdb[1];
            axyz[2]= coor_pdb[2];

        //    beg_p=rand_point; 
    //    cout<<"XXX2"<<endl;  

            int tmp3=0;
            vector<poseCoord > fin_mat;
        //    fin_mat = tmp_mat;
        //    if(beg_p==0)
        //    {
                tmp3=0;
                vector<poseCoord>().swap(fin_mat);
                for(int j=0;j<pnum;j++)
                {
                    fin_mat.push_back(pointsBx[3*j+0]);
                    fin_mat.push_back(pointsBx[3*j+1]);
                    fin_mat.push_back(pointsBx[3*j+2]);
                    tmp3=tmp3+1;
                }
        //        cout<<"XXX0"<<endl;
        /*        dis_p=0.0;
                tp_clashx= false;
                old_vwd =0.0;
                old_clash = 0.0;
                for(int jj=0;jj<fin_mat.size();jj++)
                {
                    if(fin_mat[jj].elt_=="CA")
                    {
                        for(int t=0;t<pointsBx.size();t++)
                        {
                            if(pointsBx[t].elt_=="CA")
                            {
                            //    if(abs(jj+3*beg_p-t)<1) continue;
                            //    if(abs(t)<3*end_p) continue;
                                if(abs(t)>(3*beg_p-1) && abs(t)<3*(end_p)) continue;
                                dis_p = Distance_point(fin_mat[jj].x_,pointsBx[t].x_);
                                if(dis_p < 3.70)
                                {
                                    old_clash = old_clash + 1.0/ sqrt(dis_p) ;
                //                    old_clash = old_clash + exp( 3.75-dis_p) ;
                                }
                            }
                        }
                    }
                }    
                float d1=0.0;
                float d2=0.0;
                if(randop>=0)
                {
                    d1= Distance_point(pointsBx[3*(randop-1)+1].x_,fin_mat[1].x_); // the distance between the beg_p CA atom and (beg_p-1) CA
                    d2= Distance_point(pointsBx[3*(randop-1)+1].x_,fin_mat[4].x_); // the distance between the beg_p CA atom and (beg_p-2) CA
                //    if(d2<d1) continue;
                } 
                if(randopx>=0)
                {
                    d1= Distance_point(pointsBx[3*(randopx)+1].x_,fin_mat[fin_mat.size()-2].x_); // the distance between the end_p CA atom and (end-1) CA
                    d2= Distance_point(pointsBx[3*(randopx)+1].x_,fin_mat[fin_mat.size()-5].x_); // the distance between the end_p CA atom and (end-2) CA
                //    if(d2<d1) continue;
                }                 */
                old_CC1=new_CC2;   
                old_dE = old_CC1 ;// + old_clash + 0.1*d1/d2;// + olddis_p_m;
            //    old_dE = old_clash + olddis_p_m;                      
                vector<poseCoord > fin_matx;
                fin_matx=fin_mat;
                tmp3=0;
                for(int j=0;j<fin_mat.size();j++) // not last N
                {
                    fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                    fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                    fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                    tmp3=tmp3+1;
                }
                vector<poseCoord>().swap(fin_matx);
        //        cout<<"XXX1"<<endl;
        /*        dis_p=0.0;
                new_vwd =0.0;
                new_clash = 0.0;
            //    cout<<"size: "<<fin_mat.size()<<" "<<pointsBx.size()<<endl;
                for(int jj=0;jj<fin_mat.size();jj++)
                {
                    if(fin_mat[jj].elt_=="CA")
                    {
                        for(int t=0;t<pointsBx.size();t++)
                        {
                            if(pointsBx[t].elt_=="CA")
                            {
                            //    if(abs(t)<3*end_p) continue;
                                if(abs(t)>(3*beg_p-1) && abs(t)<3*(end_p)) continue;
                                dis_p = Distance_point(fin_mat[jj].x_,pointsBx[t].x_);
                                if(dis_p < 3.70)
                                {
                                    new_clash = new_clash + 1.0/ sqrt(dis_p) ;
                 //                   cout<<"jj,t: "<<jj<<" "<<t<<" "<<dis_p<<endl;
                //                    old_clash = old_clash + exp( 3.75-dis_p) ;
                                }
                            }
                        }
                    }
                }  

        if(randop>=0)
        {
            d1= Distance_point(pointsBx[3*(randop-1)+1].x_,fin_mat[1].x_); // the distance between the beg_p CA atom and (beg_p-1) CA
            d2= Distance_point(pointsBx[3*(randop-1)+1].x_,fin_mat[4].x_); // the distance between the beg_p CA atom and (beg_p-2) CA
        //    if(d2<d1) continue;
        } 
        if(randopx>=0)
        {
            d1= Distance_point(pointsBx[3*(randopx)+1].x_,fin_mat[fin_mat.size()-2].x_); // the distance between the end_p CA atom and (end-1) CA
            d2= Distance_point(pointsBx[3*(randopx)+1].x_,fin_mat[fin_mat.size()-5].x_); // the distance between the end_p CA atom and (end-2) CA
        //    if(d2<d1) continue;
        }     */

        //    } 
        /*    else 
            {
                tmp3=0;
                vector<poseCoord>().swap(fin_mat);
                for(int j=beg_p;j<end_p;j++)
                {
                    fin_mat.push_back(pointsBx[3*j+0]);
                    fin_mat.push_back(pointsBx[3*j+1]);
                    fin_mat.push_back(pointsBx[3*j+2]);
                    tmp3=tmp3+1;
                } 
            //    cout<<"fin_mat coor0: "<<fin_mat[1].x_[0]<<" "<<fin_mat[1].x_[1]<<" "<<fin_mat[1].x_[2]<<endl;
                dis_p=0.0;
                old_vwd =0.0;
                old_clash = 0.0;
        //        cout<<"size: "<<fin_mat.size()<<" "<<pointsBx.size()<<endl;
                for(int jj=0;jj<fin_mat.size();jj++)
                {
                    if(fin_mat[jj].elt_=="CA")
                    {
                        for(int t=0;t<pointsBx.size();t++)
                        {
                            if(pointsBx[t].elt_=="CA")
                            {
                                if(abs(t)>=3*beg_p) continue;
                                dis_p = Distance_point(fin_mat[jj].x_,pointsBx[t].x_);
                                if(dis_p < 3.70)
                                {
                                    old_clash = old_clash + 3.75*1.0/ sqrt(dis_p) ;
                            //        cout<<"jj,t,dis_p: "<<jj<<" "<<t<<" "<<dis_p<<endl;
                                //    cout<<"jjd: "<<fin_mat[jj].x_[0]<<" "<<fin_mat[jj].x_[1]<<" "<<fin_mat[jj].x_[2]<<endl;
                                //    cout<<"pd: "<<pointsBx[t].x_[0]<<" "<<pointsBx[t].x_[1]<<" "<<pointsBx[t].x_[2]<<endl;
                //                    old_clash = old_clash + exp( 3.75-dis_p) ;
                                }
                            }
                        }
                    }
                }    
                old_CC1=new_CC2;   
                old_dE = old_CC1 + old_clash;//+ olddis_p_m; 
            //     old_dE =  old_clash+ olddis_p_m;                   
                vector<poseCoord > fin_matx=fin_mat;    
                tmp3=0;        
                for(int j=1;j<fin_mat.size();j++) // not first N
                {
                    fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                    fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                    fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                    tmp3=tmp3+1;
                }
                dis_p=0.0;
                new_vwd =0.0;
                new_clash = 0.0;
        //        cout<<"size: "<<fin_mat.size()<<" "<<pointsBx.size()<<endl;
                for(int jj=0;jj<fin_mat.size();jj++)
                {
                    if(fin_mat[jj].elt_=="CA")
                    {
                        for(int t=0;t<pointsBx.size();t++)
                        {
                            if(pointsBx[t].elt_=="CA")
                            {
                                if(abs(t)>=3*beg_p) continue;
                                dis_p = Distance_point(fin_mat[jj].x_,pointsBx[t].x_);
                                if(dis_p < 3.70)
                                {
                                    new_clash = new_clash + 3.75*1.0/ sqrt(dis_p) ;
                 //                   cout<<"jj,t: "<<jj<<" "<<t<<" "<<dis_p<<endl;
                //                    old_clash = old_clash + exp( 3.75-dis_p) ;
                                }
                            }
                        }
                    }
                }                                         
            }       */

        //    cout<<"fin_mat coor1: "<<fin_mat[1].x_[0]<<" "<<fin_mat[1].x_[1]<<" "<<fin_mat[1].x_[2]<<endl;  
        //    if(new_clash>10) continue;
            bond0max=-10000.0;bond1max=-10000.0;bond2max=-10000.0;
            bond0min=10000.0;bond1min=10000.0;bond2min=10000.0;            
        //    vector<vector<float> > fin_matx;
            del_ijx=vector<float>(3,0.0);
            atm_jx=vector<float>(3,0.0); 
            tmp_i=0;     
            for ( int t=0; t<theDensityMapx.density.u1()*theDensityMapx.density.u2()*theDensityMapx.density.u3(); ++t ) {
                rhoC01[t]=0.0;
                inv_rho_mask1[t]=1.0;
            }  
    //        vector<vector<float> >().swap(atm_idx); 
            atm_idx=vector<float>(3,0.0);
            for(int i=0;i<fin_mat.size();i++)  // rotate fragment from beg_p to end_p
            {
                 if(fin_mat[i].elt_ == "CA")
                {           
        //            cout<<"RRRR: "<<i;
                    vector<float> cartX1;
                    vector<float> fracX1;
                    elt_i = fin_mat[i].elt_;
                    elt_i = elt_i[0];
                    OneGaussianScattering sig_j = get_A( elt_i );
                    k = sig_j.k( theDensityMapx.effectiveB );
                    C = sig_j.C( k );
                    if ( C < 1e-6 ) continue;   

                    cartX1 = fin_mat[i].x_;           
                    MatrixTimesTransVector(theDensityMapx.c2f,cartX1,fracX1);

                    // the location of atom in grid ?
                    atm_idx[0] = pos_mod (double(fracX1[0]*theDensityMapx.grid[0] - theDensityMapx.origin[0] + 1) , (double)theDensityMapx.grid[0]);
                    atm_idx[1] = pos_mod (double(fracX1[1]*theDensityMapx.grid[1] - theDensityMapx.origin[1] + 1) , (double)theDensityMapx.grid[1]);
                    atm_idx[2] = pos_mod (double(fracX1[2]*theDensityMapx.grid[2] - theDensityMapx.origin[2] + 1) , (double)theDensityMapx.grid[2]);   
            //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
            //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
            //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
            //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
                    for(int z=1;z<=theDensityMapx.density.u3();z++)
                    {
                        atm_jx[2] = z;
                        del_ijx[2] =(atm_idx[2]-atm_jx[2])/theDensityMapx.grid[2];
                        if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                        if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                        del_ijx[0] = del_ijx[1] = 0.0;
                        vector<float> frac_tmpz;
                        MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,frac_tmpz);
                        if(square_len(frac_tmpz)> (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;                    
                        for(int y=1;y<=theDensityMapx.density.u2();y++)
                        {
                            atm_jx[1] = y;
                            del_ijx[1] = (atm_idx[1] - atm_jx[1])/theDensityMapx.grid[1];
                            // wrap-around??
                            if ( del_ijx[1] > 0.5 ) del_ijx[1]-=1.0;
                            if ( del_ijx[1] < -0.5 ) del_ijx[1]+=1.0;                        
                            del_ijx[0] = 0.0;
                            vector<float> frac_tmpy;
                            MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,frac_tmpy);
                            if(square_len(frac_tmpy)> (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;                        
                            for(int x=1;x<=theDensityMapx.density.u1();x++)
                            {
                                atm_jx[0] = x;
                                del_ijx[0] = (atm_idx[0] - atm_jx[0])/theDensityMapx.grid[0];
                                // wrap-around??
                                if ( del_ijx[0] > 0.5 ) del_ijx[0]-=1.0;
                                if ( del_ijx[0] < -0.5 ) del_ijx[0]+=1.0;                            
                                vector<float> cart_del_ij2;
                                MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,cart_del_ij2);
                                float d2 = square_len(cart_del_ij2);
                                if(d2 > (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;
                            
                                float atm = C*exp(-k*d2);
        //                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
                                float sigmoid_msk = exp( d2 - (theDensityMapx.ATOM_MASK)*(theDensityMapx.ATOM_MASK)  );
                            //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
                            //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
                                float inv_msk = 1/(1+sigmoid_msk);
                                rhoC01(x,y,z) += atm;
                                inv_rho_mask1(x,y,z) *= (1 - inv_msk);


                                if(x>bond0max) bond0max = x;
                                if(y>bond1max) bond1max = y;
                                if(z>bond2max) bond2max = z;
                                if(x<bond0min) bond0min = x;
                                if(y<bond1min) bond1min = y;
                                if(z<bond2min) bond2min = z;  
        //                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
        //                                mask2(x,y,z) = 1.0; // problem?
        //                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
        //                                    float atm = C*exp(-k*d2);
        //                                    rhoC(x,y,z) += atm;
        //                                }
        //                            }
                            }
                        }
                    } 
        //            cout<<" "<<tmp_i<<" "<<endl;           
                    tmp_i = tmp_i + 1;
                } 
            } 
    /*        ObjexxFCL::FArray3D< std::complex<double> > Finv_rho_mask2;
            fourier::fft3(inv_rho_mask1, Finv_rho_mask2);
            //int H,K,L;
            for ( int z=1; z<=(int)theDensityMap.grid[2]; ++z ) {
                int H = (z < (int)theDensityMap.grid[2]/2) ? z-1 : z-theDensityMap.grid[2] - 1;
                for ( int y=1; y<=(int)theDensityMap.grid[1]; ++y ) {
                    int K = (y < (int)theDensityMap.grid[1]/2) ? y-1 : y-theDensityMap.grid[1]-1;
                    for ( int x=1; x<=(int)theDensityMap.grid[0]; ++x ) {
                        int L = (x < (int)theDensityMap.grid[0]/2) ? x-1 : x-theDensityMap.grid[0]-1;
                        float S2c =  theDensityMap.S2(H,K,L);

                        // exp fade
                        float highres_limit=500.0;
                        float scale = exp(-S2c*(highres_limit*highres_limit));
                        Finv_rho_mask2(x,y,z) *= scale;
                    }
                }
            }
            fourier::ifft3(Finv_rho_mask2, inv_rho_mask1); */

            sumC_i2=0.0; sumO_i2=0.0; sumCO_i2=0.0; vol_i2=0.0; CC_i2=0.0;
            sumO2_i2=0.0; sumC2_i2=0.0; varC_i2=0.0; varO_i2=0.0;
            clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
            allpdbmass=0.0;
            allmrcmass=0.0;
            mass_pdb=vector<float>(3,0.0);
            mass_mrc=vector<float>(3,0.0);            
            for(int x=int(bond0min);x<=bond0max;x++)
            {
               for(int y=int(bond1min);y<=bond1max;y++)
                {
                    for(int z=int(bond2min);z<=bond2max;z++)
                    {
                        clc_x2 = rhoC01(x,y,z);
        //                clc_x2 = tmp_rhc;
                        obs_x2 = theDensityMapx.density(x,y,z);
        //                obs_x2 = tmp_den;
                        eps_x2 = 1-inv_rho_mask1(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                        sumCO_i2 += eps_x2*clc_x2*obs_x2;
                        sumO_i2  += eps_x2*obs_x2;
                        sumO2_i2 += eps_x2*obs_x2*obs_x2;
                        sumC_i2  += eps_x2*clc_x2;
                        sumC2_i2 += eps_x2*clc_x2*clc_x2;
                        vol_i2   += eps_x2;
                        mass_pdb[0]=  mass_pdb[0] + eps_x2*clc_x2*x;
                        mass_pdb[1]=  mass_pdb[1] + eps_x2*clc_x2*y;
                        mass_pdb[2]=  mass_pdb[2] + eps_x2*clc_x2*z;
                        mass_mrc[0]=  mass_mrc[0] + eps_x2*obs_x2*x;
                        mass_mrc[1]=  mass_mrc[1] + eps_x2*obs_x2*y;
                        mass_mrc[2]=  mass_mrc[2] + eps_x2*obs_x2*z;
                        allpdbmass = allpdbmass + eps_x2*clc_x2;
                        allmrcmass = allmrcmass + eps_x2*obs_x2;                        
                    }
                }
            }                       
            varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
            varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
            if ( varC_i2 == 0 || varO_i2 == 0 ) {
                CC_i2 = 0;
            } else {
                CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
            }   
            mass_pdb[0] = mass_pdb[0]/allpdbmass;
            mass_pdb[1] = mass_pdb[1]/allpdbmass;
            mass_pdb[2] = mass_pdb[2]/allpdbmass;
            mass_mrc[0] = mass_mrc[0]/allmrcmass;
            mass_mrc[1] = mass_mrc[1]/allmrcmass;
            mass_mrc[2] = mass_mrc[2]/allmrcmass;
//            newdis_p_m=Distance_point(mass_pdb,mass_mrc);             
/*            for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
    //            tmp_den=0.0;
    //            tmp_rhc =0.0;
    //            if(theDensityMap.density[x] > (2.0/3.0)*max_den) tmp_den= theDensityMap.density[x];
    //            if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
    //            if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                // fetch this point
                clc_x2 = rhoC2[x];
    //            clc_x2 = tmp_rhc;
                obs_x2 = theDensityMap.density[x];
    //            obs_x2 = tmp_den;
                eps_x2 = 1-inv_rho_mask2[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
    //            eps_x2 = 1.0;

                // SMOOTHED
                sumCO_i2 += eps_x2*clc_x2*obs_x2;
                sumO_i2  += eps_x2*obs_x2;
                sumO2_i2 += eps_x2*obs_x2*obs_x2;
                sumC_i2  += eps_x2*clc_x2;
                sumC2_i2 += eps_x2*clc_x2*clc_x2;
                vol_i2   += eps_x2;
            }
            varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
            varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
            if ( varC_i2 == 0 || varO_i2 == 0 ) {
                CC_i2 = 0;
            } else {
                CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
            }   */           
        // new clash
        //    cout<<"new_clash: "<<new_clash<<" ";           
            new_CC1=1.0-CC_i2; 
            new_dE = new_CC1 ;// + new_clash+0.1*d1/d2;// + newdis_p_m;
        //    new_dE = new_clash + newdis_p_m;
       //    cout<<"old_CC1,new_CC1:R "<<old_CC1<<" "<<new_CC1<<endl;
            cout<<"old_CC1,new_CC1:RRRRRRRRRRRRR "<<old_CC1<<" "<<new_CC1<<endl;
            dE = new_dE-old_dE;
            acc_rate=acc_rate+1;
            if(new_dE<old_dE)
            {
                tmp3=0;
                coor_pdb=vector<float>(3,0.0);
                for(int j=0;j<pnum;j++)
                {
                    pointsBx[3*j+0]=fin_mat[3*tmp3+0];
                    pointsBx[3*j+1]=fin_mat[3*tmp3+1];
                    pointsBx[3*j+2]=fin_mat[3*tmp3+2];

                    coor_pdb[0]= coor_pdb[0] + fin_mat[3*tmp3+1].x_[0];
                    coor_pdb[1]= coor_pdb[1] + fin_mat[3*tmp3+1].x_[1];
                    coor_pdb[2]= coor_pdb[2] + fin_mat[3*tmp3+1].x_[2];                    
                    tmp3=tmp3+1;
                }
                coor_pdb[0] = coor_pdb[0]/(pnum);
                coor_pdb[1] = coor_pdb[1]/(pnum);
                coor_pdb[2] = coor_pdb[2]/(pnum); 
                if(new_dE<best_E)
                {
                    best_E = new_dE;
                    vector<poseCoord>().swap(best_model);
                    best_model = fin_mat;
                }
            //    cout<<"coor: "<<coor_pdb[0]<<" "<<coor_pdb[1]<<" "<<coor_pdb[2]<<endl;                  
                new_CC2 = new_CC1;
                acc_rate=0;
                olddis_p_m = newdis_p_m;
                 
                cout<<"old_CC1,new_CC1: "<<old_dE<<" "<<new_dE<<endl;         
            }else
            {
            //    float tmpx=rand()/double(RAND_MAX);
                float tmpx=randf0and1();
                float mc_v = exp(-dE/(KT));
                if(tmpx< mc_v) // CC must be >0
                {
                    tmp3=0;
                    coor_pdb=vector<float>(3,0.0);
                    for(int j=0;j<pnum;j++)
                    {
                        pointsBx[3*j+0]=fin_mat[3*tmp3+0];
                        pointsBx[3*j+1]=fin_mat[3*tmp3+1];
                        pointsBx[3*j+2]=fin_mat[3*tmp3+2];

                        coor_pdb[0]= coor_pdb[0] + fin_mat[3*tmp3+1].x_[0];
                        coor_pdb[1]= coor_pdb[1] + fin_mat[3*tmp3+1].x_[1];
                        coor_pdb[2]= coor_pdb[2] + fin_mat[3*tmp3+1].x_[2];                         
                        tmp3=tmp3+1;
                    } 
                    coor_pdb[0] = coor_pdb[0]/(pnum);
                    coor_pdb[1] = coor_pdb[1]/(pnum);
                    coor_pdb[2] = coor_pdb[2]/(pnum); 
                //    cout<<"coor: "<<coor_pdb[0]<<" "<<coor_pdb[1]<<" "<<coor_pdb[2]<<endl;                     
                    new_CC2 = new_CC1;
                    acc_rate=0;
                    olddis_p_m = newdis_p_m;
                    cout<<"old_CC1,new_CC1:XX "<<old_dE<<" "<<new_dE<<endl;                    
                }                
            }           
        }      
    }

        old_vwd =0.0;
        old_clash = 0.0;
        old_dE =0.0;
        new_dE =0.0;
        axyz=vector<float>(3,0.0);  
        dE=0.0;  
        bond0max=-10000.0;bond1max=-10000.0;bond2max=-10000.0;
        bond0min=10000.0;bond1min=10000.0;bond2min=10000.0;            
    //    vector<vector<float> > fin_matx;
        new_CC1=0.0;
        old_CC1=0.0;
        new_CC2=0.0;
        new_dE2=0.0;
        del_ijx=vector<float>(3,0.0);
        atm_jx=vector<float>(3,0.0); 
        tmp_i=0;     
        ObjexxFCL::FArray3D< float > rhoC02;
        ObjexxFCL::FArray3D< float > inv_rho_mask2;
        rhoC02.dimension(theDensityMap.density.u1() , theDensityMap.density.u2() , theDensityMap.density.u3());
        inv_rho_mask2.dimension(theDensityMap.density.u1() , theDensityMap.density.u2() , theDensityMap.density.u3());
        for ( int t=0; t<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++t ) {
            rhoC02[t]=0.0;
            inv_rho_mask2[t]=1.0;
        }        
        vector<poseCoord > fin_maty;
        coor_pdb=vector<float>(3,0.0);
        for(int j=0;j<pnum;j++)
        {
            fin_maty.push_back(pointsBx[3*j+0]);
            fin_maty.push_back(pointsBx[3*j+1]);
            fin_maty.push_back(pointsBx[3*j+2]);
            coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
            coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
            coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];
        }
        best_model = fin_maty;
        best_E =10.0;
        coor_pdb[0] = coor_pdb[0]/(pnum);
        coor_pdb[1] = coor_pdb[1]/(pnum);
        coor_pdb[2] = coor_pdb[2]/(pnum);
        cout<<"coor: "<<coor_pdb[0]<<" "<<coor_pdb[1]<<" "<<coor_pdb[2]<<endl; 
        atm_idx=vector<float>(3,0.0);
        for(int i=0;i<fin_maty.size();i++)  // rotate fragment from beg_p to end_p
        {
             if(fin_maty[i].elt_ == "CA")
            {           
    //            cout<<"RRRR: "<<i;
                vector<float> cartX1;
                vector<float> fracX1;
                elt_i = fin_maty[i].elt_;
                elt_i = elt_i[0];
                OneGaussianScattering sig_j = get_A( elt_i );
                k = sig_j.k( theDensityMap.effectiveB );
                C = sig_j.C( k );
                if ( C < 1e-6 ) continue;   

                cartX1 = fin_maty[i].x_;           
                MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);

                // the location of atom in grid ?
                atm_idx[0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                atm_idx[1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                atm_idx[2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);   
        //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
        //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
        //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
        //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
                for(int z=1;z<=theDensityMap.density.u3();z++)
                {
                    atm_jx[2] = z;
                    del_ijx[2] =(atm_idx[2]-atm_jx[2])/theDensityMap.grid[2];
                    if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                    if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                    del_ijx[0] = del_ijx[1] = 0.0;
                    vector<float> frac_tmpz;
                    MatrixTimesTransVector(theDensityMap.f2c,del_ijx,frac_tmpz);
                    if(square_len(frac_tmpz)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                    
                    for(int y=1;y<=theDensityMap.density.u2();y++)
                    {
                        atm_jx[1] = y;
                        del_ijx[1] = (atm_idx[1] - atm_jx[1])/theDensityMap.grid[1];
                        // wrap-around??
                        if ( del_ijx[1] > 0.5 ) del_ijx[1]-=1.0;
                        if ( del_ijx[1] < -0.5 ) del_ijx[1]+=1.0;                        
                        del_ijx[0] = 0.0;
                        vector<float> frac_tmpy;
                        MatrixTimesTransVector(theDensityMap.f2c,del_ijx,frac_tmpy);
                        if(square_len(frac_tmpy)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                        
                        for(int x=1;x<=theDensityMap.density.u1();x++)
                        {
                            atm_jx[0] = x;
                            del_ijx[0] = (atm_idx[0] - atm_jx[0])/theDensityMap.grid[0];
                            // wrap-around??
                            if ( del_ijx[0] > 0.5 ) del_ijx[0]-=1.0;
                            if ( del_ijx[0] < -0.5 ) del_ijx[0]+=1.0;                            
                            vector<float> cart_del_ij2;
                            MatrixTimesTransVector(theDensityMap.f2c,del_ijx,cart_del_ij2);
                            float d2 = square_len(cart_del_ij2);
                            if(d2 > (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;
                        
                            float atm = C*exp(-k*d2);
    //                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
                            float sigmoid_msk = exp( d2 - (theDensityMap.ATOM_MASK)*(theDensityMap.ATOM_MASK)  );
                        //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
                        //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
                            float inv_msk = 1/(1+sigmoid_msk);
                            rhoC02(x,y,z) += atm;
                            inv_rho_mask2(x,y,z) *= (1 - inv_msk);


                            if(x>bond0max) bond0max = x;
                            if(y>bond1max) bond1max = y;
                            if(z>bond2max) bond2max = z;
                            if(x<bond0min) bond0min = x;
                            if(y<bond1min) bond1min = y;
                            if(z<bond2min) bond2min = z;
    //                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
    //                                mask2(x,y,z) = 1.0; // problem?
    //                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
    //                                    float atm = C*exp(-k*d2);
    //                                    rhoC(x,y,z) += atm;
    //                                }
    //                            }
                        }
                    }
                } 
    //            cout<<" "<<tmp_i<<" "<<endl;           
                tmp_i = tmp_i + 1;
            } 
        }    
          
        sumC_i2=0.0; sumO_i2=0.0; sumCO_i2=0.0; vol_i2=0.0; CC_i2=0.0;
        sumO2_i2=0.0; sumC2_i2=0.0; varC_i2=0.0; varO_i2=0.0;
        clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
        allpdbmass=0.0;
        allmrcmass=0.0;
        mass_pdb=vector<float>(3,0.0);
        mass_mrc=vector<float>(3,0.0);
        for(int x=int(bond0min);x<=bond0max;x++)
        {
           for(int y=int(bond1min);y<=bond1max;y++)
            {
                for(int z=int(bond2min);z<=bond2max;z++)
                {
                    clc_x2 = rhoC02(x,y,z);
    //                clc_x2 = tmp_rhc;
                    obs_x2 = theDensityMap.density(x,y,z);
    //                obs_x2 = tmp_den;
                    eps_x2 = 1-inv_rho_mask2(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                    sumCO_i2 += eps_x2*clc_x2*obs_x2;
                    sumO_i2  += eps_x2*obs_x2;
                    sumO2_i2 += eps_x2*obs_x2*obs_x2;
                    sumC_i2  += eps_x2*clc_x2;
                    sumC2_i2 += eps_x2*clc_x2*clc_x2;
                    vol_i2   += eps_x2;
                }
            }
        }
    //    olddis_p_m=Distance_point(mass_pdb,mass_mrc);   
        /*            for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
            //            tmp_den=0.0;
            //            tmp_rhc =0.0;
            //            if(theDensityMap.density[x] > (2.0/3.0)*max_den) tmp_den= theDensityMap.density[x];
            //            if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
            //            if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                        // fetch this point
                        clc_x2 = rhoC2[x];
            //            clc_x2 = tmp_rhc;
                        obs_x2 = theDensityMap.density[x];
            //            obs_x2 = tmp_den;
                        eps_x2 = 1-inv_rho_mask2[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
            //            eps_x2 = 1.0;

                        // SMOOTHED
                        sumCO_i2 += eps_x2*clc_x2*obs_x2;
                        sumO_i2  += eps_x2*obs_x2;
                        sumO2_i2 += eps_x2*obs_x2*obs_x2;
                        sumC_i2  += eps_x2*clc_x2;
                        sumC2_i2 += eps_x2*clc_x2*clc_x2;
                        vol_i2   += eps_x2;
                    }
                    varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
                    varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
                    if ( varC_i2 == 0 || varO_i2 == 0 ) {
                        CC_i2 = 0;
                    } else {
                        CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
                    }   */           
        varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
        varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 );
        if ( varC_i2 == 0 || varO_i2 == 0 ) {
            CC_i2 = 0;
        } else {
            CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
        }   
        new_CC2=1.0-CC_i2;   
    //    acc_rate=0;
    //    cout<<"Rend: "<<Rend[ttt]<<endl;
        KT=0.003;
        for(int iii=0;iii<200;iii++) // 200
        {
    //        cout<<"trdisx,anglex,KT: "<<trdis0n<<" "<<angle0n<<" "<<KT<<endl;
    //        float asin_theta=2.0*RandomDoubleX(0.0,1.0)-1.0;
//            float KT0 =0.1;
//            float KTn =0.001;
//            float KTx= pow(float(KT1/KT0),float(float(iii)/float(499)));
//            KT=KT0*float(KTx); 

            float asin_theta=2.0*randf0and1()-1.0;
            float acos_theta=sqrt(1.0-asin_theta*asin_theta);
            float apha=2.0*PI*randf0and1();
            float awx=acos_theta*cos(apha);
            float awy=acos_theta*sin(apha);
            float awz=asin_theta;
            // Translation Vector
            float t0=0.5;
        //    float t0=1.0;

        //    cout<<"t0: "<<t1<<" "<<t2<<" "<<t3<<endl;

            // Rotation matrix
            float anggg=10.0;
            float angle_rotategg=(2.0*randf0and1()-1.0)*anggg; // rotate angle
            float t1=(randf0and1()*2.0-1.0)*t0+0.0;
            float t2=(randf0and1()*2.0-1.0)*t0+0.0;
            float t3=(randf0and1()*2.0-1.0)*t0+0.0;  
        //    cout<<"ang,t0: "<<angle_rotategg<<" ("<<t1<<","<<t2<<","<<t3<<")"<<endl;          
            
            float asin=sin(angle_rotategg);
            float acos=cos(angle_rotategg);
            float u[3][3];
            u[0][0]=acos+awx*awx*(1.0-acos);
            u[0][1]=awx*awy*(1.0-acos)-awz*asin;
            u[0][2]=awx*awz*(1.0-acos)+awy*asin;
            u[1][0]=awx*awy*(1.0-acos)+awz*asin;
            u[1][1]=acos+awy*awy*(1.0-acos);
            u[1][2]=awy*awz*(1.0-acos)-awx*asin;
            u[2][0]=awx*awz*(1.0-acos)-awy*asin;
            u[2][1]=awy*awz*(1.0-acos)+awx*asin;
            u[2][2]=acos+awz*awz*(1.0-acos);            
            // Rotation points
    //      int rand_point = rand()%(pnumx);
    //      int rand_point = rand()%(pnum);
    //      axyz[0]=pointsx[rand_point].ca[0];
    //      axyz[1]=pointsx[rand_point].ca[1];
    //      axyz[2]=pointsx[rand_point].ca[2];
        //    int rand_point = rand()%6+83 ;
        //    int rand_point = randop ;
        //    axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
        //    axyz[1]=pointsBx[3*rand_point+1].x_[1];
        //    axyz[2]=pointsBx[3*rand_point+1].x_[2];
            axyz[0]= coor_pdb[0];  
            axyz[1]= coor_pdb[1];
            axyz[2]= coor_pdb[2];

        //    beg_p=rand_point; 
    //    cout<<"XXX2"<<endl;  

            int tmp3=0;
            vector<poseCoord > fin_mat;
        //    fin_mat = tmp_mat;
        //    if(beg_p==0)
        //    {
                tmp3=0;
                vector<poseCoord>().swap(fin_mat);
                for(int j=0;j<pnum;j++)
                {
                    fin_mat.push_back(pointsBx[3*j+0]);
                    fin_mat.push_back(pointsBx[3*j+1]);
                    fin_mat.push_back(pointsBx[3*j+2]);
                    tmp3=tmp3+1;
                }
        //        cout<<"XXX0"<<endl;
        /*        dis_p=0.0;
                tp_clashx= false;
                old_vwd =0.0;
                old_clash = 0.0;
                for(int jj=0;jj<fin_mat.size();jj++)
                {
                    if(fin_mat[jj].elt_=="CA")
                    {
                        for(int t=0;t<pointsBx.size();t++)
                        {
                            if(pointsBx[t].elt_=="CA")
                            {
                            //    if(abs(jj+3*beg_p-t)<1) continue;
                            //    if(abs(t)<3*end_p) continue;
                                if(abs(t)>(3*beg_p-1) && abs(t)<3*(end_p)) continue;
                                dis_p = Distance_point(fin_mat[jj].x_,pointsBx[t].x_);
                                if(dis_p < 3.70)
                                {
                                    old_clash = old_clash + 1.0/ sqrt(dis_p) ;
                //                    old_clash = old_clash + exp( 3.75-dis_p) ;
                                }
                            }
                        }
                    }
                }    
                float d1=0.0;
                float d2=0.0;
                if(randop>=0)
                {
                    d1= Distance_point(pointsBx[3*(randop-1)+1].x_,fin_mat[1].x_); // the distance between the beg_p CA atom and (beg_p-1) CA
                    d2= Distance_point(pointsBx[3*(randop-1)+1].x_,fin_mat[4].x_); // the distance between the beg_p CA atom and (beg_p-2) CA
                //    if(d2<d1) continue;
                } 
                if(randopx>=0)
                {
                    d1= Distance_point(pointsBx[3*(randopx)+1].x_,fin_mat[fin_mat.size()-2].x_); // the distance between the end_p CA atom and (end-1) CA
                    d2= Distance_point(pointsBx[3*(randopx)+1].x_,fin_mat[fin_mat.size()-5].x_); // the distance between the end_p CA atom and (end-2) CA
                //    if(d2<d1) continue;
                }                 */
                old_CC1=new_CC2;   
                old_dE = old_CC1 ;// + old_clash + 0.1*d1/d2;// + olddis_p_m;
            //    old_dE = old_clash + olddis_p_m;                      
                vector<poseCoord > fin_matx;
                fin_matx=fin_mat;
                tmp3=0;
                for(int j=0;j<fin_mat.size();j++) // not last N
                {
                    fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                    fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                    fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                    tmp3=tmp3+1;
                }
                vector<poseCoord>().swap(fin_matx);
        //        cout<<"XXX1"<<endl;

        //    cout<<"fin_mat coor1: "<<fin_mat[1].x_[0]<<" "<<fin_mat[1].x_[1]<<" "<<fin_mat[1].x_[2]<<endl;  
        //    if(new_clash>10) continue;
            bond0max=-10000.0;bond1max=-10000.0;bond2max=-10000.0;
            bond0min=10000.0;bond1min=10000.0;bond2min=10000.0;            
        //    vector<vector<float> > fin_matx;
            del_ijx=vector<float>(3,0.0);
            atm_jx=vector<float>(3,0.0); 
            tmp_i=0;     
            for ( int t=0; t<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++t ) {
                rhoC02[t]=0.0;
                inv_rho_mask2[t]=1.0;
            }  
//            vector<vector<float> >().swap(atm_idx); 
            atm_idx=vector<float>(3,0.0);
            for(int i=0;i<fin_mat.size();i++)  // rotate fragment from beg_p to end_p
            {
                 if(fin_mat[i].elt_ == "CA")
                {           
        //            cout<<"RRRR: "<<i;
                    vector<float> cartX1;
                    vector<float> fracX1;
                    elt_i = fin_mat[i].elt_;
                    elt_i = elt_i[0];
                    OneGaussianScattering sig_j = get_A( elt_i );
                    k = sig_j.k( theDensityMap.effectiveB );
                    C = sig_j.C( k );
                    if ( C < 1e-6 ) continue;   

                    cartX1 = fin_mat[i].x_;           
                    MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);

                    // the location of atom in grid ?
                    atm_idx[0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                    atm_idx[1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                    atm_idx[2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);   
            //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
            //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
            //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
            //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
                    for(int z=1;z<=theDensityMap.density.u3();z++)
                    {
                        atm_jx[2] = z;
                        del_ijx[2] =(atm_idx[2]-atm_jx[2])/theDensityMap.grid[2];
                        if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                        if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                        del_ijx[0] = del_ijx[1] = 0.0;
                        vector<float> frac_tmpz;
                        MatrixTimesTransVector(theDensityMap.f2c,del_ijx,frac_tmpz);
                        if(square_len(frac_tmpz)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                    
                        for(int y=1;y<=theDensityMap.density.u2();y++)
                        {
                            atm_jx[1] = y;
                            del_ijx[1] = (atm_idx[1] - atm_jx[1])/theDensityMap.grid[1];
                            // wrap-around??
                            if ( del_ijx[1] > 0.5 ) del_ijx[1]-=1.0;
                            if ( del_ijx[1] < -0.5 ) del_ijx[1]+=1.0;                        
                            del_ijx[0] = 0.0;
                            vector<float> frac_tmpy;
                            MatrixTimesTransVector(theDensityMap.f2c,del_ijx,frac_tmpy);
                            if(square_len(frac_tmpy)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                        
                            for(int x=1;x<=theDensityMap.density.u1();x++)
                            {
                                atm_jx[0] = x;
                                del_ijx[0] = (atm_idx[0] - atm_jx[0])/theDensityMap.grid[0];
                                // wrap-around??
                                if ( del_ijx[0] > 0.5 ) del_ijx[0]-=1.0;
                                if ( del_ijx[0] < -0.5 ) del_ijx[0]+=1.0;                            
                                vector<float> cart_del_ij2;
                                MatrixTimesTransVector(theDensityMap.f2c,del_ijx,cart_del_ij2);
                                float d2 = square_len(cart_del_ij2);
                                if(d2 > (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;
                            
                                float atm = C*exp(-k*d2);
        //                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
                                float sigmoid_msk = exp( d2 - (theDensityMap.ATOM_MASK)*(theDensityMap.ATOM_MASK)  );
                            //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
                            //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
                                float inv_msk = 1/(1+sigmoid_msk);
                                rhoC02(x,y,z) += atm;
                                inv_rho_mask2(x,y,z) *= (1 - inv_msk);


                                if(x>bond0max) bond0max = x;
                                if(y>bond1max) bond1max = y;
                                if(z>bond2max) bond2max = z;
                                if(x<bond0min) bond0min = x;
                                if(y<bond1min) bond1min = y;
                                if(z<bond2min) bond2min = z;  
        //                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
        //                                mask2(x,y,z) = 1.0; // problem?
        //                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
        //                                    float atm = C*exp(-k*d2);
        //                                    rhoC(x,y,z) += atm;
        //                                }
        //                            }
                            }
                        }
                    } 
        //            cout<<" "<<tmp_i<<" "<<endl;           
                    tmp_i = tmp_i + 1;
                } 
            } 

            sumC_i2=0.0; sumO_i2=0.0; sumCO_i2=0.0; vol_i2=0.0; CC_i2=0.0;
            sumO2_i2=0.0; sumC2_i2=0.0; varC_i2=0.0; varO_i2=0.0;
            clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
            allpdbmass=0.0;
            allmrcmass=0.0;
            mass_pdb=vector<float>(3,0.0);
            mass_mrc=vector<float>(3,0.0);            
            for(int x=int(bond0min);x<=bond0max;x++)
            {
               for(int y=int(bond1min);y<=bond1max;y++)
                {
                    for(int z=int(bond2min);z<=bond2max;z++)
                    {
                        clc_x2 = rhoC02(x,y,z);
        //                clc_x2 = tmp_rhc;
                        obs_x2 = theDensityMap.density(x,y,z);
        //                obs_x2 = tmp_den;
                        eps_x2 = 1-inv_rho_mask2(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                        sumCO_i2 += eps_x2*clc_x2*obs_x2;
                        sumO_i2  += eps_x2*obs_x2;
                        sumO2_i2 += eps_x2*obs_x2*obs_x2;
                        sumC_i2  += eps_x2*clc_x2;
                        sumC2_i2 += eps_x2*clc_x2*clc_x2;
                        vol_i2   += eps_x2;                    
                    }
                }
            }                       
            varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
            varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 );
            if ( varC_i2 == 0 || varO_i2 == 0 ) {
                CC_i2 = 0;
            } else {
                CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
            }   
//            newdis_p_m=Distance_point(mass_pdb,mass_mrc);             
/*            for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
    //            tmp_den=0.0;
    //            tmp_rhc =0.0;
    //            if(theDensityMap.density[x] > (2.0/3.0)*max_den) tmp_den= theDensityMap.density[x];
    //            if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
    //            if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                // fetch this point
                clc_x2 = rhoC2[x];
    //            clc_x2 = tmp_rhc;
                obs_x2 = theDensityMap.density[x];
    //            obs_x2 = tmp_den;
                eps_x2 = 1-inv_rho_mask2[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
    //            eps_x2 = 1.0;

                // SMOOTHED
                sumCO_i2 += eps_x2*clc_x2*obs_x2;
                sumO_i2  += eps_x2*obs_x2;
                sumO2_i2 += eps_x2*obs_x2*obs_x2;
                sumC_i2  += eps_x2*clc_x2;
                sumC2_i2 += eps_x2*clc_x2*clc_x2;
                vol_i2   += eps_x2;
            }
            varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
            varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
            if ( varC_i2 == 0 || varO_i2 == 0 ) {
                CC_i2 = 0;
            } else {
                CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
            }   */           
        // new clash
        //    cout<<"new_clash: "<<new_clash<<" ";           
            new_CC1=1.0-CC_i2; 
            new_dE = new_CC1 ;// + new_clash+0.1*d1/d2;// + newdis_p_m;
        //    new_dE = new_clash + newdis_p_m;
       //    cout<<"old_CC1,new_CC1:R "<<old_CC1<<" "<<new_CC1<<endl;
            cout<<"old_CC1,new_CC1:RRRRRRRRRRRRR "<<old_CC1<<" "<<new_CC1<<endl;
            dE = new_dE-old_dE;
            if(new_dE<old_dE)
            {
                tmp3=0;
                coor_pdb=vector<float>(3,0.0);
                for(int j=0;j<pnum;j++)
                {
                    pointsBx[3*j+0]=fin_mat[3*tmp3+0];
                    pointsBx[3*j+1]=fin_mat[3*tmp3+1];
                    pointsBx[3*j+2]=fin_mat[3*tmp3+2];

                    coor_pdb[0]= coor_pdb[0] + fin_mat[3*tmp3+1].x_[0];
                    coor_pdb[1]= coor_pdb[1] + fin_mat[3*tmp3+1].x_[1];
                    coor_pdb[2]= coor_pdb[2] + fin_mat[3*tmp3+1].x_[2];                    
                    tmp3=tmp3+1;
                }
                coor_pdb[0] = coor_pdb[0]/(pnum);
                coor_pdb[1] = coor_pdb[1]/(pnum);
                coor_pdb[2] = coor_pdb[2]/(pnum); 
                if(new_dE<best_E)
                {
                    best_E = new_dE;
                    vector<poseCoord>().swap(best_model);
                    best_model = fin_mat;
                }
            //    cout<<"coor: "<<coor_pdb[0]<<" "<<coor_pdb[1]<<" "<<coor_pdb[2]<<endl;                  
                new_CC2 = new_CC1;
                 
                cout<<"old_CC1,new_CC1: "<<old_dE<<" "<<new_dE<<endl;         
            }else
            {
            //    float tmpx=rand()/double(RAND_MAX);
                float tmpx=randf0and1();
                float mc_v = exp(-dE/(KT));
                if(tmpx< mc_v) // CC must be >0
                {
                    tmp3=0;
                    coor_pdb=vector<float>(3,0.0);
                    for(int j=0;j<pnum;j++)
                    {
                        pointsBx[3*j+0]=fin_mat[3*tmp3+0];
                        pointsBx[3*j+1]=fin_mat[3*tmp3+1];
                        pointsBx[3*j+2]=fin_mat[3*tmp3+2];

                        coor_pdb[0]= coor_pdb[0] + fin_mat[3*tmp3+1].x_[0];
                        coor_pdb[1]= coor_pdb[1] + fin_mat[3*tmp3+1].x_[1];
                        coor_pdb[2]= coor_pdb[2] + fin_mat[3*tmp3+1].x_[2];                         
                        tmp3=tmp3+1;
                    } 
                    coor_pdb[0] = coor_pdb[0]/(pnum);
                    coor_pdb[1] = coor_pdb[1]/(pnum);
                    coor_pdb[2] = coor_pdb[2]/(pnum); 
                //    cout<<"coor: "<<coor_pdb[0]<<" "<<coor_pdb[1]<<" "<<coor_pdb[2]<<endl;                     
                    new_CC2 = new_CC1;
                    cout<<"old_CC1,new_CC1:XX "<<old_dE<<" "<<new_dE<<endl;                    
                }                
            }           
        }          

    cout<<"XXX2"<<endl;  
        int tmp3=0;
/*        for(int j=0;j<pnum;j++)
        {
            pointsBx[3*j+0]=best_model[3*tmp3+0];
            pointsBx[3*j+1]=best_model[3*tmp3+1];
            pointsBx[3*j+2]=best_model[3*tmp3+2];              
            tmp3=tmp3+1;
        }
        cout<<"best_E: "<<best_E<<endl;    */
    } 

    cout<<"CCCCCC"<<endl;
    best_model = pointsBx;
    best_E = 100000.0;
    int NN0=0;
    int NN1=0;
    int NN2=0;
    int NN3=0;
    int NN4=0;
    int NN5=0;
    int NN6=0;
    int NN7=0;
    bool NN0TF=false;
    bool NN1TF=false;
    bool NN2TF=false;
    bool NN3TF=false;
    bool NN4TF=false;
    bool NN5TF=false;
    bool NN6TF=false;
    bool NN7TF=false;      
    // second structure rigid body movement
    bool calres = true;
    vector<float> res_CC(pnum,0.0);
    vector<sssegment> segm;
    int ssss=0;
    int rec_c = 30; // 5
    int rec_b = 20;
    int rec_a = 500;
    vector<int> num_rand;

    vector<float> res_CCk(pnum,0.0);
    vector<float> res_CC3k(pnum,0.0);
    int nhg=0;
    int numc=0;
    int nhgc=0;
    for(int j=0;j<pnum;j++)
    {
        vector<poseCoord> sing_res(1);
        sing_res[0]= pointsBx[3*j+1];
        res_CCk[j] = theDensityMap.matchposex(sing_res);
        res_CC3k[j]=res_CCk[j];
        if(res_CC3k[j]<0.05) nhg=nhg+1;
    //        cout<<"j,CC: "<<j<<" "<<res_CC[j]<<endl; 
    }
//    vector<float> res_CC3k(pnum,0.0);

/*    for(int j=0;j<pnum;j++)
    {
        for(int jp=-1;jp<2;jp++)
        {
            int jjp = j+ jp;
            if(jjp<0) jjp=0;
            if(jjp>pnum-1) jjp=pnum-1;
            res_CC3k[j] = res_CC3k[j] + res_CCk[jjp];
        }
        res_CC3k[j] = res_CC3k[j]/3.0;
        if(res_CC3k[j]<0.05) nhg=nhg+1;
    } */
    vector<distmap> dist_map;
    for(int i=0;i<pnum;i++)
    {
        if(res_CC3k[i]<0.06) continue;
        for(int j=i;j<pnum;j++)
        {
            if(res_CC3k[j]<0.06) continue;
            if(abs(i-j)<2) continue;
            float tmp_dist = Distance_point(pointsBx[3*i+1].x_,pointsBx[3*j+1].x_);
            if(tmp_dist<7.0)
    //        if(tmp_dist>4.0)
            {
                distmap tmp_map;
                tmp_map.ires = i;
                tmp_map.jres = j;
                tmp_map.dv = tmp_dist;
                dist_map.push_back(tmp_map);
            }
        }
    }
    int N_distg=dist_map.size();

    vector<distmap> dist_pre;
    string tmp_file=bindiry;
    string dist_file=tmp_file+"/"+dist_name+"/"+"dist_gauss.txt";
    ifstream myfile(dist_file);
    string tempdd;
    if(!myfile.is_open())
    {
        cout<<"Error opening distance map."<<endl;
    }else
    {
        while(getline(myfile,tempdd))
        {
            vector<string> dist_nux = string_splitxx(tempdd,' ');
            float tmp_pp=atof(dist_nux[3].c_str());
    //        if(tmp_pp<0.3) continue;
            int ig=atoi(dist_nux[0].c_str())-1;
            int jg=atoi(dist_nux[1].c_str())-1;
    //        if(res_CC3k[ig]>0.08) continue;
    //        if(res_CC3k[jg]>0.08) continue;
            distmap tmp_map;
            tmp_map.ires = atoi(dist_nux[0].c_str())-1;
            tmp_map.jres = atoi(dist_nux[1].c_str())-1;
            tmp_map.dv = atof(dist_nux[2].c_str());
            tmp_map.prob = atof(dist_nux[3].c_str());
            tmp_map.ave = atof(dist_nux[4].c_str());
            tmp_map.stdx = atof(dist_nux[5].c_str());
            dist_pre.push_back(tmp_map);
    //        cout<<"i,j,dis,prob,ave,std: "<<tmp_map.ires<<" "<<tmp_map.jres<<" "<<tmp_map.dv<<" "<<tmp_map.prob<<" "<<tmp_map.ave<<" "<<tmp_map.stdx<<endl;
        }
    }
    myfile.close();
    int N_dist=dist_pre.size();
    float w1=500.0;
    float w6=60.0;
    float w11=500.0;
    float w66=60.0;
    w1=w11*(float(pnum-nhg)/float(pnum));
    w6=w66*(float(pnum)/float(pnum-nhg)); 
//    vector<poseCoords> num_init_stru(rec_c);
//    num_init_stru.push_back(pointsBx); // initial model   
    for(int jjy=0;jjy<rec_c;jjy++)
    {
        float KT0s = 3.0 ;//0.5;
        float KT0e = 0.5 ;//0.1;
        float KT0x= pow(float(KT0e/KT0s),float(float(jjy)/float(rec_c)));
        float KT0xx = KT0s*float(KT0x);
        float KTns = 0.1;
        float KTne = 0.01;
        float KTnx= pow(float(KTne/KTns),float(float(jjy)/float(rec_c)));
        float KTnxx = KTns*float(KTnx);  
        for(int jjz=0;jjz<rec_b;jjz++)
        {          
            float KTx= pow(float(KTnxx/KT0xx),float(float(jjz)/float(rec_b)));
            float KT = KT0xx*float(KTx);   
    //        cout<<"KT: "<<KT<<endl;
                        // calculate CC between the density map and ervery residue
            vector<sssegment>().swap(segm);
            if(calres)
            {
                res_CC = vector<float>(pnum,0.0);
            //    nhg=0;
            //    vector<float> res_CC3(pnum,0.0);
                for(int j=0;j<pnum;j++)
                {
                    vector<poseCoord> sing_res(1);
                    sing_res[0]= pointsBx[3*j+1];
                    res_CC[j] = theDensityMap.matchposex(sing_res);
                    res_CC3k[j]= res_CC[j];
                //    res_CC3[j]=res_CC3k[j];
                    vector<int>().swap(num_rand);
                    if(res_CC3k[j]<0.05)
                    {
                        num_rand.push_back(j);
                    }
            //        cout<<"j,CC: "<<j<<" "<<res_CC[j]<<endl; 
                }
                
                
            /*    for(int j=0;j<pnum;j++)
                {
                    for(int jp=-1;jp<2;jp++)
                    {
                        int jjp = j+ jp;
                        if(jjp<0) jjp=0;
                        if(jjp>pnum-1) jjp=pnum-1;
                        res_CC3[j] = res_CC3[j] + res_CC[jjp];
                    }
                    res_CC3[j] = res_CC3[j]/3.0;
                    res_CC3k[j]=res_CC3[j];
                //    if(res_CC3k[j]<0.05) nhg=nhg+1;
                    vector<int>().swap(num_rand);
                    if(res_CC3[j]<0.05)
                    {
                        num_rand.push_back(j);
                    }
                } */
                for(int j=0;j<pnum;j++)
                {
                    int initx0=0;
                    int initx1=0;
                    if(res_CC3k[j]<=0.1) // 0 is the cut off value for structure in the out of density map
                    {
                        int kx0=j-1;
                        if(kx0<0) kx0=0;
                        initx0=kx0;
                        for(int i=j;i<pnum;i++)
                        {
                            if(res_CC3k[i]>0.1) 
                            {
                                float randcx=(2.0*randf0and1()-1.0);                                
                                if(randcx<0.0) break;
                        //        break;
                            }
                            j++;
                        }
                        int kx1=j+1;
                        if(kx1>pnum-1) kx1 = pnum-1;
                        initx1=kx1;
                    }
                    int dis_01 = initx1-initx0+1;
                    if(dis_01>=5) // 7 is the number which judgement how many residue should be as fragment
                    {
                        sssegment smg;
                        smg.init = initx0;
                        smg.term = initx1;
                        segm.push_back(smg);
                    }
                }
            }                 
            for(int jjx=0;jjx<rec_a;jjx++)
            {          
                vector<sssegment> segmx;
                bool TF_segmx = false;
                if(segm.size()>0)
                {
                    TF_segmx=true;
                }
                if(TF_segmx)
                {
                    int i=0;
                    int last_segmx = segm.size()-1;
                    int ia1 = segm[i].init;
                    int ia2 = segm[i].term;
                    if(ia1<3)
                    {
                        sssegment smga;   
                        smga.init = ia1;
                        smga.term = ia2; 
                        for(int j=i+1;j<=last_segmx;j++)
                        {
                            int ic1 = segm[j].init;
                            int ic2 = segm[j].term;
                            if(abs(ic1-ia2)<7)
                            {
                                smga.term = ic2;
                                ia2 = ic2;
                                i++;
                            } else
                            {                
                                break;
                            }               
                        }                  
                        segmx.push_back(smga);                         
                    } else
                    {
                        sssegment smga;   
                        smga.init = ia1;
                        smga.term = ia2;  
                        segmx.push_back(smga);                       
                    }

                    int ib1 = segm[last_segmx].init;
                    int ib2 = segm[last_segmx].term; 
                    sssegment smgb;
                    if(ib2==(pnum-1))
                    {
            //            sssegment smg;   
                        smgb.init = ib1;
                        smgb.term = ib2; 
                        for(int j=last_segmx-1;j>=0;j--)
                        {
                            int ic1 = segm[j].init;
                            int ic2 = segm[j].term;
                            if(abs(ic2-ib1)<7)
                            {
                                smgb.init = ic1;
                                ib1 = ic1;
                                last_segmx--;
                    //            i++;
                            } else
                            {                
                                break;
                            }               
                        }                  
            //            segmx.push_back(smg);                                               
                    } else
                    {
                        smgb.init = ib1;
                        smgb.term = ib2;                        
                    }
                    for(int k=i+1;k<last_segmx;k++)
                    {
                        sssegment smg;
                        smg.init = segm[k].init;
                        smg.term = segm[k].term; 
                        segmx.push_back(smg);
                    }
                    segmx.push_back(smgb);
                }
                bool TF_segmy = false;
                if(segmx.size()>0)
                {
                    TF_segmy = true;
                }                 
                float new_CC1=0.0;
                float old_CC1=0.0;
          //      float new_CC2=0.0;
         //       float KT0s = 0.01;
         //       float KT0e = 0.0001;
        //        float KT0x= pow(float(KT0e/KT0s),float(float(jjx)/float(2000)));
        //        float KT = KT0s*float(KT0x);       
        //        KT =0.001;
                ang=90.0;
                int rand_a=0,rand_b=0;
                int rand_a2=0,rand_b2=0;
                int randtx =0;
        //        if(jjy<int(0.2*float(rec_c)))
                if(jjx%5==0)
                {
                    int num_randss=num_rand.size();
                    if(num_randss>0)
                    {
                        int numrndst=randIntCustom(0,num_randss-1);
                        randtx=num_rand[numrndst];
                    } else
                    {
                        randtx = randIntCustom(0,pnum-1);
                    }
                    
                }else
                {
                    randtx = randIntCustom(0,pnum-1);
                }
                
                { 
                    int randt = randIntCustom(0,4)+3;
                    int jmx=randtx,jnx=randtx;
                    float randt_num = 2.0*randf0and1()-1.0;
                    if(randt_num>0)
                    {
                        for(int j=0;j<randt;j++)
                        {
                            jnx=randtx+j;
                            if(jnx>=(pnum-1)) break;
                        }
                        rand_a=jmx;
                        rand_b=jnx;                                                
                    }else
                    {
                        for(int j=0;j<randt;j++)
                        {
                            jmx=randtx-j;
                            if(jmx<=0) break;
                        } 
                        rand_a=jmx;
                        rand_b=jnx;                                               
                    }                
                }

                if(TF_segmy)
                {
             //       float randt_vs = 2.0*randf0and1()-1.0;
                //    if(randt_vs>=0)
                //    {
                        for(int i=0;i<segmx.size();i++)
                        {
                            int io1 = segmx[i].init;
                            int io2 = segmx[i].term;
                         //   if(rand_b<io2 && rand_a>io1) 
                         //   {
                         //       rand_a = io1;
                         //       rand_b = io2;
                         //   }                           
                            if(rand_a<io2 && rand_a>io1) rand_a = io1;
                            if(rand_b<io2 && rand_b>io1) rand_b = io2;
                        }
                //    } 
                }

                int NSS=0;
                float randt_vg = 2.0*randf0and1()-1.0;
                if(randt_vg>=0)
                {
                    for(int i=1;i<pnum;i++)
                    {
                        int rand_ai = rand_a-i;
                        if(rand_ai<0) break;
                        if(numtoss(bb[rand_ai].sst)=='H'||numtoss(bb[rand_ai].sst)=='E')
                        {
                            rand_a = rand_ai;
                        } else
                        {
                            break;
                        }
                    }
                    for(int i=1;i<pnum;i++)
                    {
                        int rand_bi = rand_b+i;
                        if(rand_bi>(pnum-1)) break;
                    //    float randt_num = 2.0*randf0and1()-1.0;
                        if(numtoss(bb[rand_bi].sst)=='H'||numtoss(bb[rand_bi].sst)=='E')
                        {
                            rand_b = rand_bi;
                        } else
                        {
                            break;
                        }               
                    }                     
                }

/*                for(int i=0;i<segm.size();i++)
                {
                    int io1 = segm[i].init;
                    int io2 = segm[i].term;
                    if(rand_a<io2 && rand_a>io1) rand_a = io1;
                    if(rand_b<io2 && rand_b>io1) rand_b = io2;
                }
                for(int i=1;i<pnum;i++)
                {
                    int rand_ai = rand_a-i;
                    if(rand_ai<0) break;
                    if(numtoss(bb[rand_a].sst)=='C') break;
                    if(numtoss(bb[rand_ai].sst)==numtoss(bb[rand_a].sst))
                    {
                        rand_a = rand_ai;
                    } else
                    {
                        break;
                    }
                }
                for(int i=1;i<pnum;i++)
                {
                    int rand_bi = rand_b+i;
                    if(rand_bi>(pnum-1)) break;
                    if(numtoss(bb[rand_b].sst)=='C') break;
                    if(numtoss(bb[rand_bi].sst)==numtoss(bb[rand_b].sst))
                    {
                        rand_b = rand_bi;
                    } else
                    {
                        break;
                    }
                } */

                rand_a2 = rand_a-2;
                rand_b2 = rand_b+2;
                if(rand_a <0) rand_a = 0;
                if(rand_b >(pnum-1)) rand_b= pnum-1;
                if(rand_a2<0) rand_a2=0;
                if(rand_b2>(pnum-1)) rand_b2= pnum-1;
            //    cout<<" "<<numtoss(bb[randtx].sst)<<endl;
            //    cout<<"rand a b: "<<rand_a<<" "<<rand_b<<endl;
                if(abs(rand_b-rand_a+1)<2) continue;
             //   cout<<"rand a b: "<<rand_a<<" "<<rand_b<<endl;
                
                vector<poseCoord > fin_maty;
                vector<poseCoord>().swap(fin_maty);
                coor_pdb=vector<float>(3,0.0);
                for(int j=0;j<pnum;j++)
                {
                    fin_maty.push_back(pointsBx[3*j+0]);
                    fin_maty.push_back(pointsBx[3*j+1]);
                    fin_maty.push_back(pointsBx[3*j+2]);
                }     
                nhgc=0;   
                for(int j=rand_a;j<=rand_b;j++)
                {
                    coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
                    coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
                    coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];  
                    if(res_CC3k[j]<0.05) nhgc=nhgc+1;           
                }
                coor_pdb[0] = coor_pdb[0]/abs(rand_b-rand_a+1);
                coor_pdb[1] = coor_pdb[1]/abs(rand_b-rand_a+1);
                coor_pdb[2] = coor_pdb[2]/abs(rand_b-rand_a+1);  
                numc=rand_b-rand_a+1;
                w1=w11*(float(numc-nhgc)/float(numc));
                w6=w66*(float(numc)/float(numc-nhgc)); 
                if(w1==0)
                {
                    w1=0.5*w11;
                    w6=3.0*w66;
                } 
                del_ijx=vector<float>(3,0.0);
                atm_jx=vector<float>(3,0.0); 
                tmp_i=0;     
                for ( int t=0; t<pDenxyz; ++t ) {
                    rhoC0x[t]=0.0;
                    inv_rho_mask0x[t]=1.0;
                }
                atm_idx=vector<float>(3,0.0);
                bond0max=-1000.0;bond1max=-1000.0;bond2max=-1000.0;
                bond0min=1000.0;bond1min=1000.0;bond2min=1000.0;      
                for(int i=3*rand_a2;i<=3*rand_b2+2;i++)  // rotate fragment from beg_p to end_p
                {
                     if(fin_maty[i].elt_ == "CA")
                    {           
                        vector<float> cartX1;
                        vector<float> fracX1;
                        elt_i = fin_maty[i].elt_;
                        elt_i = elt_i[0];
                        OneGaussianScattering sig_j = get_A( elt_i );
                        k = sig_j.k( theDensityMapx.effectiveB );
                        C = sig_j.C( k );
                        if ( C < 1e-6 ) continue;   

                        cartX1 = fin_maty[i].x_;           
                        MatrixTimesTransVector(theDensityMapx.c2f,cartX1,fracX1);

                        // the location of atom in grid ?
                        atm_idx[0] = pos_mod (double(fracX1[0]*pgrid[0] - porigin[0] + 1) , (double)pgrid[0]);
                        atm_idx[1] = pos_mod (double(fracX1[1]*pgrid[1] - porigin[1] + 1) , (double)pgrid[1]);
                        atm_idx[2] = pos_mod (double(fracX1[2]*pgrid[2] - porigin[2] + 1) , (double)pgrid[2]);   
                //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
                //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
                //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
                //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
                        for(int z=1;z<=pDenz;z++)
                        {
                            atm_jx[2] = z;
                            del_ijx[2] =(atm_idx[2]-atm_jx[2])/pgrid[2];
                            if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                            if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                            del_ijx[0] = del_ijx[1] = 0.0;
                            vector<float> frac_tmpz;
                            MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,frac_tmpz);
                            if(square_len(frac_tmpz)> (ppadding+patom_m)*(ppadding+patom_m) ) continue; 
                            if(z<bond2min) bond2min = z; 
                            if(z>bond2max) bond2max = z;                  
                            for(int y=1;y<=pDeny;y++)
                            {
                                atm_jx[1] = y;
                                del_ijx[1] = (atm_idx[1] - atm_jx[1])/pgrid[1];
                                // wrap-around??
                                if ( del_ijx[1] > 0.5 ) del_ijx[1]-=1.0;
                                if ( del_ijx[1] < -0.5 ) del_ijx[1]+=1.0;                        
                                del_ijx[0] = 0.0;
                                vector<float> frac_tmpy;
                                MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,frac_tmpy);
                                if(square_len(frac_tmpy)> (ppadding+patom_m)*(ppadding+patom_m) ) continue; 
                                if(y>bond1max) bond1max = y;      
                                if(y<bond1min) bond1min = y;                 
                                for(int x=1;x<=pDenx;x++)
                                {
                                    atm_jx[0] = x;
                                    del_ijx[0] = (atm_idx[0] - atm_jx[0])/pgrid[0];
                                    // wrap-around??
                                    if ( del_ijx[0] > 0.5 ) del_ijx[0]-=1.0;
                                    if ( del_ijx[0] < -0.5 ) del_ijx[0]+=1.0;                            
                                    vector<float> cart_del_ij2;
                                    MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,cart_del_ij2);
                                    float d2 = square_len(cart_del_ij2);
                                    if(d2 > (ppadding+patom_m)*(ppadding+patom_m) ) continue;
                                
                                    float atm = C*exp(-k*d2);
            //                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
                                    float sigmoid_msk = exp( d2 - (patom_m)*(patom_m)  );
                                //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
                                //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
                                    float inv_msk = 1/(1+sigmoid_msk);
                                    rhoC0x(x,y,z) += atm;
                                    inv_rho_mask0x(x,y,z) *= (1 - inv_msk);


                                    if(x>bond0max) bond0max = x;
                                    
                                    
                                    if(x<bond0min) bond0min = x;
                                    
                                    
            //                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
            //                                mask2(x,y,z) = 1.0; // problem?
            //                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
            //                                    float atm = C*exp(-k*d2);
            //                                    rhoC(x,y,z) += atm;
            //                                }
            //                            }
                                }
                            }
                        } 
            //            cout<<" "<<tmp_i<<" "<<endl;           
            //            tmp_i = tmp_i + 1;
                    } 
                } 
                float sumC_i2=0.0, sumO_i2=0.0, sumCO_i2=0.0, vol_i2=0.0, CC_i2=0.0;
                float sumO2_i2=0.0, sumC2_i2=0.0, varC_i2=0.0, varO_i2=0.0;
                float clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
        //        allpdbmass=0.0;
        //        allmrcmass=0.0;
        //        mass_pdb=vector<float>(3,0.0);
        //        mass_mrc=vector<float>(3,0.0);
                for(int x=int(bond0min);x<=bond0max;x++)
                {
                   for(int y=int(bond1min);y<=bond1max;y++)
                    {
                        for(int z=int(bond2min);z<=bond2max;z++)
                        {
                            clc_x2 = rhoC0x(x,y,z);
            //                clc_x2 = tmp_rhc;
                            obs_x2 = theDensityMapx.density(x,y,z);
            //                obs_x2 = tmp_den;
                            eps_x2 = 1-inv_rho_mask0x(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                            sumCO_i2 += eps_x2*clc_x2*obs_x2;
                            sumO_i2  += eps_x2*obs_x2;
                            sumO2_i2 += eps_x2*obs_x2*obs_x2;
                            sumC_i2  += eps_x2*clc_x2;
                            sumC2_i2 += eps_x2*clc_x2*clc_x2;
                            vol_i2   += eps_x2;
            //                mass_pdb[0]=  mass_pdb[0] + eps_x2*clc_x2*x;
            //                mass_pdb[1]=  mass_pdb[1] + eps_x2*clc_x2*y;
            //                mass_pdb[2]=  mass_pdb[2] + eps_x2*clc_x2*z;
            //                mass_mrc[0]=  mass_mrc[0] + eps_x2*obs_x2*x;
            //                mass_mrc[1]=  mass_mrc[1] + eps_x2*obs_x2*y;
            //                mass_mrc[2]=  mass_mrc[2] + eps_x2*obs_x2*z;
            //                allpdbmass = allpdbmass + eps_x2*clc_x2;
            //                allmrcmass = allmrcmass + eps_x2*obs_x2;
                        }
                    }
                }
            //    mass_pdb[0] = mass_pdb[0]/allpdbmass;
            //    mass_pdb[1] = mass_pdb[1]/allpdbmass;
            //    mass_pdb[2] = mass_pdb[2]/allpdbmass;
            //    mass_mrc[0] = mass_mrc[0]/allmrcmass;
            //    mass_mrc[1] = mass_mrc[1]/allmrcmass;
            //    mass_mrc[2] = mass_mrc[2]/allmrcmass;
            //    olddis_p_m=Distance_point(mass_pdb,mass_mrc);     
                varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
                varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 );
                if ( varC_i2 == 0 || varO_i2 == 0 ) {
                    CC_i2 = 0;
                } else {
                    CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
                }   
                old_CC1=1.0-CC_i2;
                
                float dis_p=0.0;
                old_vwd =0.0;
                old_clash = 0.0;
                float clash_s=1.0;
                for(int jj=3*rand_a2;jj<=3*rand_b2+2;jj++)
                {
                    if(fin_maty[jj].elt_=="CA")
                    {
                        string VR1 = fin_maty[jj].elt_;
                        for(int t=0;t<fin_maty.size();t++)
                        {
                            if(fin_maty[t].elt_=="CA")
                            {
                            //    if(t>=3*rand_a && t<=3*rand_b+2) continue;
                                if(abs(t-jj)<2) continue;
                            //    if(abs(t)<3*end_p) continue;
                                string VR2 = fin_maty[t].elt_ ;
                                dis_p = Distance_point(fin_maty[jj].x_,fin_maty[t].x_);
                                old_clash = old_clash + GetVdwEgCG(VR1,VR2,dis_p);
                        //        if(dis_p < 3.70)
                        //        {
                        //            old_clash = old_clash + 1.0/ sqrt(dis_p + 0.0001) ;
                        //            old_clash = old_clash + clash_s;
                                //    cout<<"jj t: "<<jj<<" "<<t<<" "<<old_clash<<endl;
                                //    cout<<"coor0: "<<fin_maty[jj].x_[0]<<" "<<fin_maty[jj].x_[1]<<" "<<fin_maty[jj].x_[2]<<endl;
                                //    cout<<"coor1: "<<pointsBx[t].x_[0]<<" "<<pointsBx[t].x_[1]<<" "<<pointsBx[t].x_[2]<<endl;
                //                    old_clash = old_clash + exp( 3.75-dis_p) ;
                        //        }
                            }
                        }
                    }
                }  

                float old_Ehbond=0.0;
        //        if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E')
                {
            /*        for(int j=0;j<pnum;j++)
                    {
                        bb[j].indn = 3*j+0;
                        bb[j].indca = 3*j+1;
                        bb[j].indc = 3*j+2;
                        string resng = pointsrenm[j];
                        char resny = resng[0];
                        bb[j].resid = aminoid(resny);                             
                    }     
                //    calcsse2(bb, pnum, pointsBx);  
                    calcssennhoc(bb,pnum,pointsBx);           */   
                    int numseq=pnum;
                    vector<point3f> decstr(numseq);                
                    tmp_tt=0;
                    for(int j=0;j<pnum;j++)
                    {
                        decstr[tmp_tt].ss2 = numtoss(bb[j].sst);// predicted secondary structure
                        decstr[tmp_tt].ssm = numtoss(bb[j].sst);
                        decstr[tmp_tt].x = pointsBx[3*tmp_tt+1].x_[0];
                        decstr[tmp_tt].y = pointsBx[3*tmp_tt+1].x_[1];
                        decstr[tmp_tt].z = pointsBx[3*tmp_tt+1].x_[2];
                        decstr[tmp_tt].ptn.x = pointsBx[3*tmp_tt+0].x_[0];
                        decstr[tmp_tt].ptn.y = pointsBx[3*tmp_tt+0].x_[1];
                        decstr[tmp_tt].ptn.z = pointsBx[3*tmp_tt+0].x_[2];
                        decstr[tmp_tt].ptc.x = pointsBx[3*tmp_tt+2].x_[0];
                        decstr[tmp_tt].ptc.y = pointsBx[3*tmp_tt+2].x_[1];
                        decstr[tmp_tt].ptc.z = pointsBx[3*tmp_tt+2].x_[2];
                        tmp_tt=tmp_tt+1;
                    } 
                //    old_Ehbond=energyhbondcanc(decstr,numseq);
                    old_Ehbond=energyhbondnhoc2(decstr,numseq);
                    vector<point3f>().swap(decstr);  
                }                 
            /*    if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E')
                {
                    for(int j=0;j<pnum;j++)
                    {
                        bb[j].indn = 3*j+0;
                        bb[j].indca = 3*j+1;
                        bb[j].indc = 3*j+2;
                        string resng = pointsrenm[j];
                        char resny = resng[0];
                        bb[j].resid = aminoid(resny);     
                    }     
                //    calcsse2(bb, pnum, pointsBx);
                    calcssennhoc(bb,pnum,pointsBx);                
                    int numseq=pnum;
                    vector<point3f> decstr(numseq);                
                    tmp_tt=0;
                    for(int j=0;j<pnum;j++)
                    {
                        decstr[tmp_tt].ss2 = numtoss(bb[j].sst);
                        decstr[tmp_tt].x = pointsBx[3*tmp_tt+1].x_[0];
                        decstr[tmp_tt].y = pointsBx[3*tmp_tt+1].x_[1];
                        decstr[tmp_tt].z = pointsBx[3*tmp_tt+1].x_[2];
                        decstr[tmp_tt].ptn.x = pointsBx[3*tmp_tt+0].x_[0];
                        decstr[tmp_tt].ptn.y = pointsBx[3*tmp_tt+0].x_[1];
                        decstr[tmp_tt].ptn.z = pointsBx[3*tmp_tt+0].x_[2];
                        decstr[tmp_tt].ptc.x = pointsBx[3*tmp_tt+2].x_[0];
                        decstr[tmp_tt].ptc.y = pointsBx[3*tmp_tt+2].x_[1];
                        decstr[tmp_tt].ptc.z = pointsBx[3*tmp_tt+2].x_[2];
                        tmp_tt=tmp_tt+1;
                    } 
            //        old_Ehbond=energyhbondcanc(decstr,numseq);
                    old_Ehbond=energyhbondnhoc2(decstr,numseq);
                    vector<point3f>().swap(decstr);  
                }  */
        /*        vector<point3f> decstrbondang(pnum);
                for(int j=0;j<pnum;j++)
                {
                    decstrbondang[j].ss2=numtoss(bb[j].sst);
                    decstrbondang[j].x=pointsBx[3*j+1].x_[0];
                    decstrbondang[j].y=pointsBx[3*j+1].x_[1];
                    decstrbondang[j].z=pointsBx[3*j+1].x_[2];
                    decstrbondang[j].ptn.x=pointsBx[3*j+0].x_[0];
                    decstrbondang[j].ptn.y=pointsBx[3*j+0].x_[1];
                    decstrbondang[j].ptn.z=pointsBx[3*j+0].x_[2];
                    decstrbondang[j].ptc.x=pointsBx[3*j+2].x_[0];
                    decstrbondang[j].ptc.y=pointsBx[3*j+2].x_[1];
                    decstrbondang[j].ptc.z=pointsBx[3*j+2].x_[2];
                    string resng = pointsrenm[j];
                    char resny = resng[0];
                    decstrbondang[j].iaa = aminoid(resny);
    //                decstr[i].iaa=pp.aminoid(decstr[i].aaa);                
    //                tmp_tt=tmp_tt+1;
                }     
                for(int j=1;j<pnum;j++)
                {        
                    if(j==0)
                    {
                        decstrbondang[j].ang[0]=0.0;
                        decstrbondang[j].ang[1]=0.0;
                        decstrbondang[j].ang[2]=0.0;
                    }                   
                    point3d p12 = setv(decstrbondang[j-1].x-decstrbondang[j-1].ptc.x,decstrbondang[j-1].y-decstrbondang[j-1].ptc.y,decstrbondang[j-1].z-decstrbondang[j-1].ptc.z);
                    point3d p23 = setv(decstrbondang[j].ptn.x-decstrbondang[j-1].ptc.x,decstrbondang[j].ptn.y-decstrbondang[j-1].ptc.y,decstrbondang[j].ptn.z-decstrbondang[j-1].ptc.z);
                    decstrbondang[j].ang[0]= angv(p12,p23)*degrad;
                    p12 = setv(decstrbondang[j-1].ptc.x-decstrbondang[j].ptn.x,decstrbondang[j-1].ptc.y-decstrbondang[j].ptn.y,decstrbondang[j-1].ptc.z-decstrbondang[j].ptn.z);
                    p23 = setv(decstrbondang[j].x-decstrbondang[j].ptn.x,decstrbondang[j].y-decstrbondang[j].ptn.y,decstrbondang[j].z-decstrbondang[j].ptn.z);
                    decstrbondang[j].ang[1]= angv(p12,p23)*degrad;
                    p12 = setv(decstrbondang[j].ptn.x-decstrbondang[j].x,decstrbondang[j].ptn.y-decstrbondang[j].y,decstrbondang[j].ptn.z-decstrbondang[j].z);
                    p23 = setv(decstrbondang[j].ptc.x-decstrbondang[j].x,decstrbondang[j].ptc.y-decstrbondang[j].y,decstrbondang[j].ptc.z-decstrbondang[j].z);
                    decstrbondang[j].ang[2]= angv(p12,p23)*degrad;
    //                tmp_tt=tmp_tt+1;
                }
        //        cout<<"FFF"<<endl; 
                float old_bondangenergy =0.0;
                Eenergybondangle(decstrbondang,pnum,old_bondangenergy);
                cout<< "old_bondangenergy: "<< old_bondangenergy/pnum<<endl;                 */                
            //    old_dE = old_CC1 + old_clash;// + olddis_p_m;        
            //    old_dE = 200.0*old_CC1 + 10.0*old_clash/abs(rand_b2-rand_a2+1) + 200.0*(old_Ehbond/abs(pnum))/(1.0-(old_Ehbond/abs(pnum)));// + old_bondangenergy/pnum;
                vector<point3f> decstrbondlen(pnum);
                for(int j=0;j<pnum;j++)
                {
                    decstrbondlen[j].ss2=numtoss(bb[j].sst);
                    decstrbondlen[j].ssm=numtoss(bb[j].sst);
                    decstrbondlen[j].stype=numtoss(bb[j].sst);
                    decstrbondlen[j].x=pointsBx[3*j+1].x_[0];
                    decstrbondlen[j].y=pointsBx[3*j+1].x_[1];
                    decstrbondlen[j].z=pointsBx[3*j+1].x_[2];
                    decstrbondlen[j].ptn.x=pointsBx[3*j+0].x_[0];
                    decstrbondlen[j].ptn.y=pointsBx[3*j+0].x_[1];
                    decstrbondlen[j].ptn.z=pointsBx[3*j+0].x_[2];
                    decstrbondlen[j].ptc.x=pointsBx[3*j+2].x_[0];
                    decstrbondlen[j].ptc.y=pointsBx[3*j+2].x_[1];
                    decstrbondlen[j].ptc.z=pointsBx[3*j+2].x_[2];
                    string resng = pointsrenm[j];
                    char resny = resng[0];
                    decstrbondlen[j].iaa = aminoid(resny);
    //                decstr[i].iaa=pp.aminoid(decstr[i].aaa);                
    //                tmp_tt=tmp_tt+1;
                } 
                str2tor(decstrbondlen,pnum,3);
            //    str2torp(decstrbondlen,pnum,0,pnum-1); 
                // bondangle energy
        /*        float old_fang=0.0;
                float old_Ebondang=0.0;
            //    cout<<"jjk: "<<endl;
                old_Ebondang =(float) Eenergybondangle(decstrbondlen,pnum,old_fang);
            //    cout<<"jjk1: "<<endl;
                old_Ebondang = 0.30*old_Ebondang + old_fang;         */
                // bondlength energy                      
                float old_fene=0.0;
                float old_Ebondlen=0.0;
                old_Ebondlen = (float) energybondlength(decstrbondlen,pnum,old_fene);
                old_Ebondlen = 0.5*old_Ebondlen + 20.0*old_fene;    

                float E_dist = 0.0;
                for(int iu=0;iu<N_dist;iu++)
                {
                    int ik = dist_pre[iu].ires;
                    int jk = dist_pre[iu].jres;
                    float tmp_x = dist_pre[iu].dv;
                    float tmp_y = Distance_point(fin_maty[3*ik+1].x_,fin_maty[3*jk+1].x_);
            //        float tmp_z=sqrt((tmp_x-tmp_y)*(tmp_x-tmp_y))-dist_pre[iu].ave;
                    float tmp_z=tmp_y-dist_pre[iu].ave;
                    E_dist = E_dist + (log(tmp_z*tmp_z+1))/sqrt(dist_pre[iu].stdx);

                }
                E_dist=E_dist/float(pnum*(rand_b2-rand_a2+1));

                float Eres0=0.0;
                for(int j=3*rand_a2;j<=3*rand_b2+2;j++)
                {
                    for(int js=j-3;js<=j+3;js++)
                    {
                        if(js==j) continue;
                        if(js<0) continue;
                        if(js>(pnum-1)) continue;
                        float dx0 = Distance_point(fin_maty[j].x_,fin_maty[js].x_);
                        float dx1 = Distance_point(pointsBx0[j].x_,pointsBx0[js].x_);
                        Eres0 = Eres0 + (dx0-dx1)*(dx0-dx1);
                    }
                }                        
                Eres0 = Eres0/float(rand_b2-rand_a2+1);

                float E_con = 0.0;
                for(int iu=0;iu<N_distg;iu++)
                {
                    int ik = dist_map[iu].ires;
                    int jk = dist_map[iu].jres;
                    float tmp_x = dist_map[iu].dv;
                    float tmp_y = Distance_point(fin_maty[3*ik+1].x_,fin_maty[3*jk+1].x_);
                    E_con = E_con + sqrt((tmp_x-tmp_y)*(tmp_x-tmp_y));
                }
                E_con=E_con/float((rand_b2-rand_a2+1));

                // torsion angle
                float old_tor_fene=0.0;
                float old_tor_E=0.0;
                old_tor_E = energyrama(decstrbondlen,pnum,old_tor_fene,rlogduke,ramaduke); 
                old_tor_E = 4.00*old_tor_E + old_tor_fene;
                old_dE = w1*old_CC1 + old_clash + old_Ebondlen + old_tor_E +1.0*E_con+10.0*Eres0+ w6*E_dist + 0.5*old_Ehbond;// + old_Ebondang/1000.0;// + 200.0*(old_Ehbond/abs(pnum));
            //    old_dE = old_CC1 + old_clash;

             /*   int rotx_point=0; // rotate point
                if(rand_a==0||rand_a==randop||rand_b==(randop-1)||rand_b==(pnum-1))
                {
                    if(rand_a==0||rand_a==randop)
                    {
                        rotx_point = rand_b;
                    }else
                    {
                        rotx_point = rand_a;
                    }

                    float asin_theta=2.0*RandomDoubleX(0.0,1.0)-1.0;
                    float acos_theta=sqrt(1.0-asin_theta*asin_theta);
                    float apha=2.0*PI*RandomDoubleX(0.0,1.0);
                    float awx=acos_theta*cos(apha);
                    float awy=acos_theta*sin(apha);
                    float awz=asin_theta;
                    // Translation Vector
                    float t0=0.1;
                    float t1=(RandomDoubleX(0.0,1.0)*2.0-1.0)*t0+0.0;
                    float t2=(RandomDoubleX(0.0,1.0)*2.0-1.0)*t0+0.0;
                    float t3=(RandomDoubleX(0.0,1.0)*2.0-1.0)*t0+0.0;
                //    cout<<"t0: "<<t1<<" "<<t2<<" "<<t3<<endl;

                    // Rotation matrix
                    float anggg=1.0;
                    float angle_rotategg=(2.0*RandomDoubleX(0.0,1.0)-1.0)*anggg; // rotate angle
                    float asin=sin(angle_rotategg);
                    float acos=cos(angle_rotategg);
                    float u[3][3];
                    u[0][0]=acos+awx*awx*(1.0-acos);
                    u[0][1]=awx*awy*(1.0-acos)-awz*asin;
                    u[0][2]=awx*awz*(1.0-acos)+awy*asin;
                    u[1][0]=awx*awy*(1.0-acos)+awz*asin;
                    u[1][1]=acos+awy*awy*(1.0-acos);
                    u[1][2]=awy*awz*(1.0-acos)-awx*asin;
                    u[2][0]=awx*awz*(1.0-acos)-awy*asin;
                    u[2][1]=awy*awz*(1.0-acos)+awx*asin;
                    u[2][2]=acos+awz*awz*(1.0-acos);       

                    int rand_point = rotx_point;
                    axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
                    axyz[1]=pointsBx[3*rand_point+1].x_[1];
                    axyz[2]=pointsBx[3*rand_point+1].x_[2];     

                    int tmp3=0;
                    vector<poseCoord > fin_mat;
                    vector<poseCoord>().swap(fin_mat);
                    fin_mat = fin_maty;
                    if(rand_a==0||rand_a==randop)
                    {
                        for(int j=0;j<fin_mat.size()-1;j++) // not last N
                        {
                            fin_mat[tmp3].x_[0]=t1+axyz[0]+(fin_maty[j].x_[0]-axyz[0])*u[0][0]+(fin_maty[j].x_[1]-axyz[1])*u[0][1]+(fin_maty[j].x_[2]-axyz[2])*u[0][2];
                            fin_mat[tmp3].x_[1]=t2+axyz[1]+(fin_maty[j].x_[0]-axyz[0])*u[1][0]+(fin_maty[j].x_[1]-axyz[1])*u[1][1]+(fin_maty[j].x_[2]-axyz[2])*u[1][2];
                            fin_mat[tmp3].x_[2]=t3+axyz[2]+(fin_maty[j].x_[0]-axyz[0])*u[2][0]+(fin_maty[j].x_[1]-axyz[1])*u[2][1]+(fin_maty[j].x_[2]-axyz[2])*u[2][2];
                            tmp3=tmp3+1;
                        }                
                    } else
                    {
                        for(int j=1;j<fin_mat.size();j++) // not first N
                        {
                            fin_mat[tmp3].x_[0]=t1+axyz[0]+(fin_maty[j].x_[0]-axyz[0])*u[0][0]+(fin_maty[j].x_[1]-axyz[1])*u[0][1]+(fin_maty[j].x_[2]-axyz[2])*u[0][2];
                            fin_mat[tmp3].x_[1]=t2+axyz[1]+(fin_maty[j].x_[0]-axyz[0])*u[1][0]+(fin_maty[j].x_[1]-axyz[1])*u[1][1]+(fin_maty[j].x_[2]-axyz[2])*u[1][2];
                            fin_mat[tmp3].x_[2]=t3+axyz[2]+(fin_maty[j].x_[0]-axyz[0])*u[2][0]+(fin_maty[j].x_[1]-axyz[1])*u[2][1]+(fin_maty[j].x_[2]-axyz[2])*u[2][2];
                            tmp3=tmp3+1;
                        }                
                    }           
                    vector<poseCoord>().swap(fin_maty);
                    fin_maty=fin_mat;                

                }
                else  */
                if(rand_b==(pnum-1)||rand_a==(0))
                {
                    int rand_movement = randIntCustom(0,1000000)%6;
                    if(rand_movement == 0)
                    {
                        float asin_theta=2.0*randf0and1()-1.0;
                        float acos_theta=sqrt(1.0-asin_theta*asin_theta);
                        float apha=2.0*PI*randf0and1();
                        float awx=acos_theta*cos(apha);
                        float awy=acos_theta*sin(apha);
                        float awz=asin_theta;
                        // Translation Vector
                        float t0=0.3;
                        float t1=(randf0and1()*2.0-1.0)*t0+0.0;
                        float t2=(randf0and1()*2.0-1.0)*t0+0.0;
                        float t3=(randf0and1()*2.0-1.0)*t0+0.0;
                    //    cout<<"t0: "<<t1<<" "<<t2<<" "<<t3<<endl;

                        // Rotation matrix
                        float anggg=180.0;
                        float angle_rotategg=(2.0*randf0and1()-1.0)*anggg; // rotate angle
                        float asin=sin(angle_rotategg*(PI/180.0));
                        float acos=cos(angle_rotategg*(PI/180.0));
                        float u[3][3];
                        u[0][0]=acos+awx*awx*(1.0-acos);
                        u[0][1]=awx*awy*(1.0-acos)-awz*asin;
                        u[0][2]=awx*awz*(1.0-acos)+awy*asin;
                        u[1][0]=awx*awy*(1.0-acos)+awz*asin;
                        u[1][1]=acos+awy*awy*(1.0-acos);
                        u[1][2]=awy*awz*(1.0-acos)-awx*asin;
                        u[2][0]=awx*awz*(1.0-acos)-awy*asin;
                        u[2][1]=awy*awz*(1.0-acos)+awx*asin;
                        u[2][2]=acos+awz*awz*(1.0-acos);
                        int rand_point = rand_a;
                        int tmp3=0;
                        vector<poseCoord > fin_mat;
                    //    fin_mat = tmp_mat;
                        tmp3=0;
                        vector<poseCoord>().swap(fin_mat);
                        fin_mat = fin_maty;
                    //    old_dE = old_clash + olddis_p_m;                      
                        vector<poseCoord > fin_matx=fin_mat;
                        tmp3=0;
                        if(rand_b==(pnum-1))
                        {
                            rand_point = rand_a;
                            axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
                            axyz[1]=pointsBx[3*rand_point+1].x_[1];
                            axyz[2]=pointsBx[3*rand_point+1].x_[2];
                            for(int j=3*rand_a+1;j<=3*rand_b+2;j++) // not last N
                            {
                                fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                                fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                                fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];                           
                                tmp3=tmp3+1;
                            } 

                        } else
                        {
                            rand_point = rand_b;
                            axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
                            axyz[1]=pointsBx[3*rand_point+1].x_[1];
                            axyz[2]=pointsBx[3*rand_point+1].x_[2];
                            for(int j=3*rand_a;j<=3*rand_b+1;j++) // not last N
                            {
                                fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                                fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                                fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                                tmp3=tmp3+1;
                            }
                                                     
                        }
                        vector<poseCoord>().swap(fin_maty);
                        fin_maty = fin_mat;                    
                        vector<poseCoord >().swap(fin_matx);
                    //    vector<poseCoord >().swap(fin_maty);
                    //    fin_maty = fin_mat;
                        vector<poseCoord >().swap(fin_mat);   
                        NN2TF = true ;                                                        
                    }
                    if(rand_movement == 1)
                    {
                        float asin_theta=2.0*randf0and1()-1.0;
                        float acos_theta=sqrt(1.0-asin_theta*asin_theta);
                        float apha=2.0*PI*randf0and1();
                        float awx=acos_theta*cos(apha);
                        float awy=acos_theta*sin(apha);
                        float awz=asin_theta;
                        // Translation Vector
                        float t0=0.3;
                        float t1=(randf0and1()*2.0-1.0)*t0+0.0;
                        float t2=(randf0and1()*2.0-1.0)*t0+0.0;
                        float t3=(randf0and1()*2.0-1.0)*t0+0.0;
                    //    cout<<"t0: "<<t1<<" "<<t2<<" "<<t3<<endl;

                        // Rotation matrix
                        float anggg=90.0;
                        float angle_rotategg=(2.0*randf0and1()-1.0)*anggg; // rotate angle
                        float asin=sin(angle_rotategg*(PI/180.0));
                        float acos=cos(angle_rotategg*(PI/180.0));
                        float u[3][3];
                        u[0][0]=acos+awx*awx*(1.0-acos);
                        u[0][1]=awx*awy*(1.0-acos)-awz*asin;
                        u[0][2]=awx*awz*(1.0-acos)+awy*asin;
                        u[1][0]=awx*awy*(1.0-acos)+awz*asin;
                        u[1][1]=acos+awy*awy*(1.0-acos);
                        u[1][2]=awy*awz*(1.0-acos)-awx*asin;
                        u[2][0]=awx*awz*(1.0-acos)-awy*asin;
                        u[2][1]=awy*awz*(1.0-acos)+awx*asin;
                        u[2][2]=acos+awz*awz*(1.0-acos);            
                        // Rotation points
                //      int rand_point = rand()%(pnumx);
                //      int rand_point = rand()%(pnum);
                //      axyz[0]=pointsx[rand_point].ca[0];
                //      axyz[1]=pointsx[rand_point].ca[1];
                //      axyz[2]=pointsx[rand_point].ca[2];
                    //    int rand_point = rand()%6+83 ;
                    //    int rand_point = randop ;
                    //    axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
                    //    axyz[1]=pointsBx[3*rand_point+1].x_[1];
                    //    axyz[2]=pointsBx[3*rand_point+1].x_[2];
                        coor_pdb=vector<float>(3,0.0);
                        for(int j=rand_a;j<=rand_b;j++)
                        {
                            coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
                            coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
                            coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];             
                        }
                        axyz[0]= coor_pdb[0]/float(rand_b-rand_a+1);  
                        axyz[1]= coor_pdb[1]/float(rand_b-rand_a+1);
                        axyz[2]= coor_pdb[2]/float(rand_b-rand_a+1);
                //        cout<<"cent: "<< randtx<<" "<< axyz[0]<<" "<<axyz[1]<<" "<< axyz[2]<<endl; 

                    //    beg_p=rand_point;   

                        int tmp3=0;
                        vector<poseCoord > fin_mat;
                    //    fin_mat = tmp_mat;
                        tmp3=0;
                        vector<poseCoord>().swap(fin_mat);
                        fin_mat = fin_maty;
                    //    old_dE = old_clash + olddis_p_m;                      
                        vector<poseCoord > fin_matx=fin_mat;
                        tmp3=0;
                        for(int j=3*rand_a;j<=3*rand_b+2;j++) // not last N
                        {
                            fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                            fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                            fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                            tmp3=tmp3+1;
                        }
                        vector<poseCoord >().swap(fin_matx);
                        int numseq=pnum;
                    //    if((rand_b-rand_a-3)<2) continue;
                    //    if(rand_a==0) continue;
                        vector<point3f> decstr(numseq);

                        for(int j=0;j<pnum;j++)
                        {
                            decstr[j].ss2=numtoss(bb[j].sst);
                            decstr[j].ssm=numtoss(bb[j].sst);
                            decstr[j].stype=numtoss(bb[j].sst);
                            decstr[j].x=fin_mat[3*j+1].x_[0];
                            decstr[j].y=fin_mat[3*j+1].x_[1];
                            decstr[j].z=fin_mat[3*j+1].x_[2];
                            decstr[j].ptn.x=fin_mat[3*j+0].x_[0];
                            decstr[j].ptn.y=fin_mat[3*j+0].x_[1];
                            decstr[j].ptn.z=fin_mat[3*j+0].x_[2];
                            decstr[j].ptc.x=fin_mat[3*j+2].x_[0];
                            decstr[j].ptc.y=fin_mat[3*j+2].x_[1];
                            decstr[j].ptc.z=fin_mat[3*j+2].x_[2];
            //                tmp_tt=tmp_tt+1;
                        }
                        
                    //    if(rand_b+2>pnum-1) continue;
                        if(rand_b==(pnum-1))
                        {
                            if((rand_a2)<1||(rand_a+3)>(pnum-1)) continue;
                            mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+3);
                        } else
                        {
                            if((rand_b-2)<1||rand_b2>(pnum-1)) continue;
                            mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b2);
                        }
                        tmp_tt=0;
                        int tmp_ttx=0;
                        for(int j=0;j<pnum;j++)
                        {
                            fin_maty[3*j+1].x_[0]=decstr[j].x;
                            fin_maty[3*j+1].x_[1]=decstr[j].y;
                            fin_maty[3*j+1].x_[2]=decstr[j].z;
                            fin_maty[3*j+0].x_[0]=decstr[j].ptn.x;
                            fin_maty[3*j+0].x_[1]=decstr[j].ptn.y;
                            fin_maty[3*j+0].x_[2]=decstr[j].ptn.z;
                            fin_maty[3*j+2].x_[0]=decstr[j].ptc.x;
                            fin_maty[3*j+2].x_[1]=decstr[j].ptc.y;
                            fin_maty[3*j+2].x_[2]=decstr[j].ptc.z;
                            tmp_tt=tmp_tt+1;
                        }
                        vector<point3f>().swap(decstr);   
                        NN1TF = true;                                      
                    }
                    if(rand_movement==2)
                    {
            //            if((rand_b+1)>(pnum-1)) continue;
                        vector<float> sftv(3,0.0);
                        sftv[0] = (2.0*randf0and1()-1.0);
                        sftv[1] = (2.0*randf0and1()-1.0);
                        sftv[2] = (2.0*randf0and1()-1.0);
                    //    cout<<"sftv: "<<sftv[0]<<" "<<sftv[1]<<" "<<sftv[2]<<endl;
                        for(int j=rand_a;j<=rand_b;j++)
                        {
                            fin_maty[3*j+0].x_[0] = fin_maty[3*(j)+0].x_[0] + sftv[0];
                            fin_maty[3*j+0].x_[1] = fin_maty[3*(j)+0].x_[1] + sftv[1];
                            fin_maty[3*j+0].x_[2] = fin_maty[3*(j)+0].x_[2] + sftv[2];
                            fin_maty[3*j+1].x_[0] = fin_maty[3*(j)+1].x_[0] + sftv[0];
                            fin_maty[3*j+1].x_[1] = fin_maty[3*(j)+1].x_[1] + sftv[1];
                            fin_maty[3*j+1].x_[2] = fin_maty[3*(j)+1].x_[2] + sftv[2];
                            fin_maty[3*j+2].x_[0] = fin_maty[3*(j)+2].x_[0] + sftv[0];
                            fin_maty[3*j+2].x_[1] = fin_maty[3*(j)+2].x_[1] + sftv[1];
                            fin_maty[3*j+2].x_[2] = fin_maty[3*(j)+2].x_[2] + sftv[2];
                        }
                        int numseq=pnum;
                    //    if((rand_b-rand_a-3)<2) continue;
                    //    if(rand_a==0) continue;
                        vector<point3f> decstr(numseq);

                        for(int j=0;j<pnum;j++)
                        {
                            decstr[j].ss2=numtoss(bb[j].sst);
                            decstr[j].ssm=numtoss(bb[j].sst);
                            decstr[j].stype=numtoss(bb[j].sst);
                            decstr[j].x=fin_maty[3*j+1].x_[0];
                            decstr[j].y=fin_maty[3*j+1].x_[1];
                            decstr[j].z=fin_maty[3*j+1].x_[2];
                            decstr[j].ptn.x=fin_maty[3*j+0].x_[0];
                            decstr[j].ptn.y=fin_maty[3*j+0].x_[1];
                            decstr[j].ptn.z=fin_maty[3*j+0].x_[2];
                            decstr[j].ptc.x=fin_maty[3*j+2].x_[0];
                            decstr[j].ptc.y=fin_maty[3*j+2].x_[1];
                            decstr[j].ptc.z=fin_maty[3*j+2].x_[2];
            //                tmp_tt=tmp_tt+1;
                        }
                        if(rand_b==(pnum-1))
                        {
                            if((rand_a2)<1||(rand_a+3)>(pnum-1)) continue;
                            mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+3);
                        } else
                        {
                            if((rand_b2)>(pnum-1)||(rand_b-2)<1) continue;
                            mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b2); 
                        }
                        for(int j=0;j<pnum;j++)
                        {
                            fin_maty[3*j+1].x_[0]=decstr[j].x;
                            fin_maty[3*j+1].x_[1]=decstr[j].y;
                            fin_maty[3*j+1].x_[2]=decstr[j].z;
                            fin_maty[3*j+0].x_[0]=decstr[j].ptn.x;
                            fin_maty[3*j+0].x_[1]=decstr[j].ptn.y;
                            fin_maty[3*j+0].x_[2]=decstr[j].ptn.z;
                            fin_maty[3*j+2].x_[0]=decstr[j].ptc.x;
                            fin_maty[3*j+2].x_[1]=decstr[j].ptc.y;
                            fin_maty[3*j+2].x_[2]=decstr[j].ptc.z;
                        }     
                        vector<point3f>().swap(decstr);  
                        NN3TF = true;              
                    } 
                    if(rand_movement==3)
                    {
            //            if((rand_b+1)>(pnum-1)) continue;
                        float rand_mov0=(randf0and1()*2.0-1.0);
                        if(rand_mov0>=0.0)
                        {
                            vector<float> sftv(3,0.0);
                            sftv[0] = fin_maty[3*rand_b+1].x_[0] -fin_maty[3*rand_a+1].x_[0];
                            sftv[1] = fin_maty[3*rand_b+1].x_[1] -fin_maty[3*rand_a+1].x_[1];
                            sftv[2] = fin_maty[3*rand_b+1].x_[2] -fin_maty[3*rand_a+1].x_[2];
                            sftv[0] = sftv[0]/float(3*rand_b-3*rand_a+1);
                            sftv[1] = sftv[1]/float(3*rand_b-3*rand_a+1);
                            sftv[2] = sftv[2]/float(3*rand_b-3*rand_a+1);
                    //        cout<<"sftv: "<<sftv[0]<<" "<<sftv[1]<<" "<<sftv[2]<<endl;
//                            sftv[0] = (2.0*randf0and1()-1.0);
//                            sftv[1] = (2.0*randf0and1()-1.0);
//                            sftv[2] = (2.0*randf0and1()-1.0);
                            for(int j=rand_a;j<=rand_b;j++)
                            {
                                fin_maty[3*j+0].x_[0] = fin_maty[3*(j)+0].x_[0] + sftv[0];
                                fin_maty[3*j+0].x_[1] = fin_maty[3*(j)+0].x_[1] + sftv[1];
                                fin_maty[3*j+0].x_[2] = fin_maty[3*(j)+0].x_[2] + sftv[2];
                                fin_maty[3*j+1].x_[0] = fin_maty[3*(j)+1].x_[0] + sftv[0];
                                fin_maty[3*j+1].x_[1] = fin_maty[3*(j)+1].x_[1] + sftv[1];
                                fin_maty[3*j+1].x_[2] = fin_maty[3*(j)+1].x_[2] + sftv[2];
                                fin_maty[3*j+2].x_[0] = fin_maty[3*(j)+2].x_[0] + sftv[0];
                                fin_maty[3*j+2].x_[1] = fin_maty[3*(j)+2].x_[1] + sftv[1];
                                fin_maty[3*j+2].x_[2] = fin_maty[3*(j)+2].x_[2] + sftv[2];
                            }
                            int numseq=pnum;
                        //    if((rand_b-rand_a-3)<2) continue;
                        //    if(rand_a==0) continue;
                            vector<point3f> decstr(numseq);

                            for(int j=0;j<pnum;j++)
                            {
                                decstr[j].ss2=numtoss(bb[j].sst);
                                decstr[j].ssm=numtoss(bb[j].sst);
                                decstr[j].stype=numtoss(bb[j].sst);
                                decstr[j].x=fin_maty[3*j+1].x_[0];
                                decstr[j].y=fin_maty[3*j+1].x_[1];
                                decstr[j].z=fin_maty[3*j+1].x_[2];
                                decstr[j].ptn.x=fin_maty[3*j+0].x_[0];
                                decstr[j].ptn.y=fin_maty[3*j+0].x_[1];
                                decstr[j].ptn.z=fin_maty[3*j+0].x_[2];
                                decstr[j].ptc.x=fin_maty[3*j+2].x_[0];
                                decstr[j].ptc.y=fin_maty[3*j+2].x_[1];
                                decstr[j].ptc.z=fin_maty[3*j+2].x_[2];
                //                tmp_tt=tmp_tt+1;
                            }
                        //    if(rand_b>(pnum-3)) 
                        //    {
                        //        int rand_a_a = rand_a-2;
                        //        int rand_a_b = rand_a+2;
                        //        if(rand_a+2>pnum-1) rand_a_b = pnum-1;
                        //        if((rand_a-2)<1) rand_a_a = 1;
                        //        mcfragsweepLMP2(decstr,numseq,rand_a_a,rand_a_b);
                        //        mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+1);
                        //    } else
                        //    {
                        //        int rand_b_a = rand_b-2;
                        //        int rand_b_b = rand_b+2;
                        //        if((rand_b-2)<1) rand_b_a=1;
                        //        if(rand_b+2>pnum-1) rand_b_b = pnum-1;
                        //        mcfragsweepLMP2(decstr,numseq,rand_b_a,rand_b_b);
                        //        mcfragsweepLMP2(decstr,numseq,rand_b-1,rand_b2); 
                        //    }
                            if(rand_b==(pnum-1))
                            {
                                if((rand_a2)<1||(rand_a+3)>(pnum-1)) continue;
                                mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+3);
                            } else
                            {
                                if((rand_b-2)<1||rand_b2>(pnum-1)) continue;
                                mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b2);
                            }                            
                            for(int j=0;j<pnum;j++)
                            {
                                fin_maty[3*j+1].x_[0]=decstr[j].x;
                                fin_maty[3*j+1].x_[1]=decstr[j].y;
                                fin_maty[3*j+1].x_[2]=decstr[j].z;
                                fin_maty[3*j+0].x_[0]=decstr[j].ptn.x;
                                fin_maty[3*j+0].x_[1]=decstr[j].ptn.y;
                                fin_maty[3*j+0].x_[2]=decstr[j].ptn.z;
                                fin_maty[3*j+2].x_[0]=decstr[j].ptc.x;
                                fin_maty[3*j+2].x_[1]=decstr[j].ptc.y;
                                fin_maty[3*j+2].x_[2]=decstr[j].ptc.z;
                            }     
                            vector<point3f>().swap(decstr);
                        } else
                        {
                            vector<float> sftv(3,0.0);
                            sftv[0] = -fin_maty[3*rand_b+1].x_[0] + fin_maty[3*rand_a+1].x_[0];
                            sftv[1] = -fin_maty[3*rand_b+1].x_[1] + fin_maty[3*rand_a+1].x_[1];
                            sftv[2] = -fin_maty[3*rand_b+1].x_[2] + fin_maty[3*rand_a+1].x_[2];
                            sftv[0] = sftv[0]/float(3*rand_b-3*rand_a+1);
                            sftv[1] = sftv[1]/float(3*rand_b-3*rand_a+1);
                            sftv[2] = sftv[2]/float(3*rand_b-3*rand_a+1);
                    //        cout<<"sftv: "<<sftv[0]<<" "<<sftv[1]<<" "<<sftv[2]<<endl;
//                            sftv[0] = (2.0*randf0and1()-1.0);
//                            sftv[1] = (2.0*randf0and1()-1.0);
//                            sftv[2] = (2.0*randf0and1()-1.0);
                            for(int j=rand_a;j<=rand_b;j++)
                            {
                                fin_maty[3*j+0].x_[0] = fin_maty[3*(j)+0].x_[0] + sftv[0];
                                fin_maty[3*j+0].x_[1] = fin_maty[3*(j)+0].x_[1] + sftv[1];
                                fin_maty[3*j+0].x_[2] = fin_maty[3*(j)+0].x_[2] + sftv[2];
                                fin_maty[3*j+1].x_[0] = fin_maty[3*(j)+1].x_[0] + sftv[0];
                                fin_maty[3*j+1].x_[1] = fin_maty[3*(j)+1].x_[1] + sftv[1];
                                fin_maty[3*j+1].x_[2] = fin_maty[3*(j)+1].x_[2] + sftv[2];
                                fin_maty[3*j+2].x_[0] = fin_maty[3*(j)+2].x_[0] + sftv[0];
                                fin_maty[3*j+2].x_[1] = fin_maty[3*(j)+2].x_[1] + sftv[1];
                                fin_maty[3*j+2].x_[2] = fin_maty[3*(j)+2].x_[2] + sftv[2];
                            }
                            int numseq=pnum;
                        //    if((rand_b-rand_a-3)<2) continue;
                        //    if(rand_a==0) continue;
                            vector<point3f> decstr(numseq);

                            for(int j=0;j<pnum;j++)
                            {
                                decstr[j].ss2=numtoss(bb[j].sst);
                                decstr[j].ssm=numtoss(bb[j].sst);
                                decstr[j].stype=numtoss(bb[j].sst);
                                decstr[j].x=fin_maty[3*j+1].x_[0];
                                decstr[j].y=fin_maty[3*j+1].x_[1];
                                decstr[j].z=fin_maty[3*j+1].x_[2];
                                decstr[j].ptn.x=fin_maty[3*j+0].x_[0];
                                decstr[j].ptn.y=fin_maty[3*j+0].x_[1];
                                decstr[j].ptn.z=fin_maty[3*j+0].x_[2];
                                decstr[j].ptc.x=fin_maty[3*j+2].x_[0];
                                decstr[j].ptc.y=fin_maty[3*j+2].x_[1];
                                decstr[j].ptc.z=fin_maty[3*j+2].x_[2];
                //                tmp_tt=tmp_tt+1;
                            }
                        //    if(rand_b>(pnum-2)) 
                        //    {
                        //        int rand_a_a = rand_a-2;
                        //        int rand_a_b = rand_a+2;
                        //        if(rand_a+2>pnum-1) rand_a_b = pnum-1;
                        //        if((rand_a-2)<1) rand_a_a = 1;
                        //        mcfragsweepLMP2(decstr,numseq,rand_a_a,rand_a_b);
                        //        mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+1);
                        //    } else
                        //    {
                        //        int rand_b_a = rand_b-2;
                        //        int rand_b_b = rand_b+2;
                        //        if((rand_b-2)<1) rand_b_a=1;
                        //        if(rand_b+2>pnum-1) rand_b_b = pnum-1;
                        //        mcfragsweepLMP2(decstr,numseq,rand_b_a,rand_b_b);
                        //        mcfragsweepLMP2(decstr,numseq,rand_b-1,rand_b2); 
                        //    }
                            if(rand_b==(pnum-1))
                            {
                                if((rand_a2)<1||(rand_a+3)>(pnum-1)) continue;
                                mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+3);
                            } else
                            {
                                if((rand_b-2)<1||rand_b2>(pnum-1)) continue;
                                mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b2);
                            }                            
                            for(int j=0;j<pnum;j++)
                            {
                                fin_maty[3*j+1].x_[0]=decstr[j].x;
                                fin_maty[3*j+1].x_[1]=decstr[j].y;
                                fin_maty[3*j+1].x_[2]=decstr[j].z;
                                fin_maty[3*j+0].x_[0]=decstr[j].ptn.x;
                                fin_maty[3*j+0].x_[1]=decstr[j].ptn.y;
                                fin_maty[3*j+0].x_[2]=decstr[j].ptn.z;
                                fin_maty[3*j+2].x_[0]=decstr[j].ptc.x;
                                fin_maty[3*j+2].x_[1]=decstr[j].ptc.y;
                                fin_maty[3*j+2].x_[2]=decstr[j].ptc.z;
                            }     
                            vector<point3f>().swap(decstr);
                        }
                        NN4TF = true;
                    }
                    if(rand_movement==4)  // random change the coordinate from atom rand_a to atom atom rand_b
                    {
            //            int numseq=abs(rand_b-rand_a+1);
                        int numseq=pnum;
                //        if((rand_b-rand_a-3)<2) continue;
                        if(rand_a==0) rand_a=1;
                        vector<point3f> decstr(numseq);

                        for(int j=0;j<pnum;j++)
                        {
                            decstr[j].ss2=numtoss(bb[j].sst);
                            decstr[j].ssm=numtoss(bb[j].sst);
                            decstr[j].stype=numtoss(bb[j].sst);
                            decstr[j].x=pointsBx[3*j+1].x_[0];
                            decstr[j].y=pointsBx[3*j+1].x_[1];
                            decstr[j].z=pointsBx[3*j+1].x_[2];
                            decstr[j].ptn.x=pointsBx[3*j+0].x_[0];
                            decstr[j].ptn.y=pointsBx[3*j+0].x_[1];
                            decstr[j].ptn.z=pointsBx[3*j+0].x_[2];
                            decstr[j].ptc.x=pointsBx[3*j+2].x_[0];
                            decstr[j].ptc.y=pointsBx[3*j+2].x_[1];
                            decstr[j].ptc.z=pointsBx[3*j+2].x_[2];
            //                tmp_tt=tmp_tt+1;
                        }                               

                    //    mcfragsweepLMP2(decstr,numseq,rand_a,rand_b);
                        mcmovementLMP(decstr,numseq,rand_a,rand_b);
                        tmp_tt=0;
                        int tmp_ttx=0;
                        for(int j=0;j<pnum;j++)
                        {
                            fin_maty[3*tmp_tt+1].x_[0]=decstr[j].x;
                            fin_maty[3*tmp_tt+1].x_[1]=decstr[j].y;
                            fin_maty[3*tmp_tt+1].x_[2]=decstr[j].z;
                            fin_maty[3*tmp_tt+0].x_[0]=decstr[j].ptn.x;
                            fin_maty[3*tmp_tt+0].x_[1]=decstr[j].ptn.y;
                            fin_maty[3*tmp_tt+0].x_[2]=decstr[j].ptn.z;
                            fin_maty[3*tmp_tt+2].x_[0]=decstr[j].ptc.x;
                            fin_maty[3*tmp_tt+2].x_[1]=decstr[j].ptc.y;
                            fin_maty[3*tmp_tt+2].x_[2]=decstr[j].ptc.z;
                            tmp_tt=tmp_tt+1;
                        }   
                        vector<point3f>().swap(decstr);  
                        NN6TF=true;           
                    } 
                    if(rand_movement==5) // roatate around axis 
                    {
                        angle_rotate = (2.0*randf0and1()-1.0)*90.0; // rotate angle
                        GroupRotationpid(fin_maty[3*rand_a+1].x_,fin_maty[3*rand_b+2-1].x_,angle_rotate,fin_maty,3*rand_a+1,3*rand_b+2-1);  
                        NN0TF =true;                     
                    }                                                            
                } else
                {
                        int rand_movement = randIntCustom(0,1000000)%6;
                   //     cout<<"sssssssssssssssss: "<<rand_movement<<endl;
                        if(rand_movement==0) // roatate around axis 
                        {
                            angle_rotate = (2.0*randf0and1()-1.0)*90.0; // rotate angle
                            GroupRotationpid(fin_maty[3*rand_a+1].x_,fin_maty[3*rand_b+2-1].x_,angle_rotate,fin_maty,3*rand_a+1,3*rand_b+2-1);  
                            NN0TF =true;                     
                        }
                        if(rand_movement==1)  // random change the coordinate from atom rand_a to atom atom rand_b
                        {
                //            int numseq=abs(rand_b-rand_a+1);
                            int numseq=pnum;
                    //        if((rand_b-rand_a-3)<2) continue;
                    //        if(rand_a==0) rand_a=1;
                            vector<point3f> decstr(numseq);

                            for(int j=0;j<pnum;j++)
                            {
                                decstr[j].ss2=numtoss(bb[j].sst);
                                decstr[j].ssm=numtoss(bb[j].sst);
                                decstr[j].stype=numtoss(bb[j].sst);
                                decstr[j].x=pointsBx[3*j+1].x_[0];
                                decstr[j].y=pointsBx[3*j+1].x_[1];
                                decstr[j].z=pointsBx[3*j+1].x_[2];
                                decstr[j].ptn.x=pointsBx[3*j+0].x_[0];
                                decstr[j].ptn.y=pointsBx[3*j+0].x_[1];
                                decstr[j].ptn.z=pointsBx[3*j+0].x_[2];
                                decstr[j].ptc.x=pointsBx[3*j+2].x_[0];
                                decstr[j].ptc.y=pointsBx[3*j+2].x_[1];
                                decstr[j].ptc.z=pointsBx[3*j+2].x_[2];
                //                tmp_tt=tmp_tt+1;
                            }                               

                        //    mcfragsweepLMP2(decstr,numseq,rand_a,rand_b);
                            mcmovementLMP(decstr,numseq,rand_a,rand_b);
                            tmp_tt=0;
                            int tmp_ttx=0;
                            for(int j=0;j<pnum;j++)
                            {
                                fin_maty[3*tmp_tt+1].x_[0]=decstr[j].x;
                                fin_maty[3*tmp_tt+1].x_[1]=decstr[j].y;
                                fin_maty[3*tmp_tt+1].x_[2]=decstr[j].z;
                                fin_maty[3*tmp_tt+0].x_[0]=decstr[j].ptn.x;
                                fin_maty[3*tmp_tt+0].x_[1]=decstr[j].ptn.y;
                                fin_maty[3*tmp_tt+0].x_[2]=decstr[j].ptn.z;
                                fin_maty[3*tmp_tt+2].x_[0]=decstr[j].ptc.x;
                                fin_maty[3*tmp_tt+2].x_[1]=decstr[j].ptc.y;
                                fin_maty[3*tmp_tt+2].x_[2]=decstr[j].ptc.z;
                                tmp_tt=tmp_tt+1;
                            }   
                            vector<point3f>().swap(decstr);  
                            NN6TF=true;           
                        }
                        if(rand_movement == 2)
                        {
                    //        if(numtoss(bb[randtx].sst)=='E') continue;
                    //        cout<<"trdisx,anglex,KT: "<<trdis0n<<" "<<angle0n<<" "<<KT<<endl;
                            float asin_theta=2.0*randf0and1()-1.0;
                            float acos_theta=sqrt(1.0-asin_theta*asin_theta);
                            float apha=2.0*PI*randf0and1();
                            float awx=acos_theta*cos(apha);
                            float awy=acos_theta*sin(apha);
                            float awz=asin_theta;
                            // Translation Vector
                            float t0=0.3;
                            float t1=(randf0and1()*2.0-1.0)*t0+0.0;
                            float t2=(randf0and1()*2.0-1.0)*t0+0.0;
                            float t3=(randf0and1()*2.0-1.0)*t0+0.0;
                        //    cout<<"t0: "<<t1<<" "<<t2<<" "<<t3<<endl;

                            // Rotation matrix
                            float anggg=90.0;
                            float angle_rotategg=(2.0*randf0and1()-1.0)*anggg; // rotate angle
                            float asin=sin(angle_rotategg*(PI/180.0));
                            float acos=cos(angle_rotategg*(PI/180.0));
                            float u[3][3];
                            u[0][0]=acos+awx*awx*(1.0-acos);
                            u[0][1]=awx*awy*(1.0-acos)-awz*asin;
                            u[0][2]=awx*awz*(1.0-acos)+awy*asin;
                            u[1][0]=awx*awy*(1.0-acos)+awz*asin;
                            u[1][1]=acos+awy*awy*(1.0-acos);
                            u[1][2]=awy*awz*(1.0-acos)-awx*asin;
                            u[2][0]=awx*awz*(1.0-acos)-awy*asin;
                            u[2][1]=awy*awz*(1.0-acos)+awx*asin;
                            u[2][2]=acos+awz*awz*(1.0-acos);            
                            // Rotation points
                    //      int rand_point = rand()%(pnumx);
                    //      int rand_point = rand()%(pnum);
                    //      axyz[0]=pointsx[rand_point].ca[0];
                    //      axyz[1]=pointsx[rand_point].ca[1];
                    //      axyz[2]=pointsx[rand_point].ca[2];
                        //    int rand_point = rand()%6+83 ;
                        //    int rand_point = randop ;
                        //    axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
                        //    axyz[1]=pointsBx[3*rand_point+1].x_[1];
                        //    axyz[2]=pointsBx[3*rand_point+1].x_[2];
                            coor_pdb=vector<float>(3,0.0);
                            for(int j=rand_a;j<=rand_b;j++)
                            {
                                coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
                                coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
                                coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];
                            }                        
                            axyz[0]= coor_pdb[0]/float(rand_b-rand_a+1);  
                            axyz[1]= coor_pdb[1]/float(rand_b-rand_a+1);
                            axyz[2]= coor_pdb[2]/float(rand_b-rand_a+1);
                    //        cout<<"cent: "<< randtx<<" "<< axyz[0]<<" "<<axyz[1]<<" "<< axyz[2]<<endl; 

                        //    beg_p=rand_point;   

                            int tmp3=0;
                            vector<poseCoord > fin_mat;
                        //    fin_mat = tmp_mat;
                            tmp3=0;
                            vector<poseCoord>().swap(fin_mat);
                            fin_mat = fin_maty;
                    /*        for(int j=rand_a+2;j<=rand_b-2;j++)
                            {
                                fin_mat.push_back(pointsBx[3*j+0]);
                                fin_mat.push_back(pointsBx[3*j+1]);
                                fin_mat.push_back(pointsBx[3*j+2]);               
                                tmp3=tmp3+1;
                            }        */
                        //    old_dE = old_clash + olddis_p_m;                      
                            vector<poseCoord > fin_matx=fin_mat;
                            tmp3=0;
                            for(int j=3*rand_a;j<=3*rand_b+2;j++) // not last N
                            {
                                fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                                fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                                fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                                tmp3=tmp3+1;
                            }
                            vector<poseCoord >().swap(fin_matx);
                            int numseq=pnum;
                        //    if((rand_b-rand_a-3)<2) continue;
                        //    if(rand_a==0) continue;
                            vector<point3f> decstr(numseq);

                            for(int j=0;j<pnum;j++)
                            {
                                decstr[j].ss2=numtoss(bb[j].sst);
                                decstr[j].ssm=numtoss(bb[j].sst);
                                decstr[j].stype=numtoss(bb[j].sst);
                                decstr[j].x=fin_mat[3*j+1].x_[0];
                                decstr[j].y=fin_mat[3*j+1].x_[1];
                                decstr[j].z=fin_mat[3*j+1].x_[2];
                                decstr[j].ptn.x=fin_mat[3*j+0].x_[0];
                                decstr[j].ptn.y=fin_mat[3*j+0].x_[1];
                                decstr[j].ptn.z=fin_mat[3*j+0].x_[2];
                                decstr[j].ptc.x=fin_mat[3*j+2].x_[0];
                                decstr[j].ptc.y=fin_mat[3*j+2].x_[1];
                                decstr[j].ptc.z=fin_mat[3*j+2].x_[2];
                //                tmp_tt=tmp_tt+1;
                            }
                            if(rand_a2<1) continue;
                            if(rand_a+2>(pnum-1)) continue;
                            if(rand_b-2<1) continue;
                            if(rand_b2>pnum-1) continue;
                            mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+2);
                            mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b2);
                            tmp_tt=0;
                            int tmp_ttx=0;
                            for(int j=0;j<pnum;j++)
                            {
                                fin_maty[3*j+1].x_[0]=decstr[j].x;
                                fin_maty[3*j+1].x_[1]=decstr[j].y;
                                fin_maty[3*j+1].x_[2]=decstr[j].z;
                                fin_maty[3*j+0].x_[0]=decstr[j].ptn.x;
                                fin_maty[3*j+0].x_[1]=decstr[j].ptn.y;
                                fin_maty[3*j+0].x_[2]=decstr[j].ptn.z;
                                fin_maty[3*j+2].x_[0]=decstr[j].ptc.x;
                                fin_maty[3*j+2].x_[1]=decstr[j].ptc.y;
                                fin_maty[3*j+2].x_[2]=decstr[j].ptc.z;
                                tmp_tt=tmp_tt+1;
                            }
                            vector<point3f>().swap(decstr); 
                            NN1TF=true;
                        } 
                        if(rand_movement == 3)
                        {
                            // Rotation axis 
                            float asin_theta=2.0*randf0and1()-1.0;
                            float acos_theta=sqrt(1.0-asin_theta*asin_theta);
                            float apha=2.0*PI*randf0and1();
                            float awx=acos_theta*cos(apha);
                            float awy=acos_theta*sin(apha);
                            float awz=asin_theta;
                            // Translation Vector
                            float t0=0.3;
                            float t1=(randf0and1()*2.0-1.0)*t0+0.0;
                            float t2=(randf0and1()*2.0-1.0)*t0+0.0;
                            float t3=(randf0and1()*2.0-1.0)*t0+0.0;

                            // Rotation matrix
                            float anggg=90.0;
                            float angle_rotategg=(2.0*randf0and1()-1.0)*anggg; // rotate angle
                            float asin=sin(angle_rotategg*(PI/180.0));
                            float acos=cos(angle_rotategg*(PI/180.0));
                            float u[3][3];
                            u[0][0]=acos+awx*awx*(1.0-acos);
                            u[0][1]=awx*awy*(1.0-acos)-awz*asin;
                            u[0][2]=awx*awz*(1.0-acos)+awy*asin;
                            u[1][0]=awx*awy*(1.0-acos)+awz*asin;
                            u[1][1]=acos+awy*awy*(1.0-acos);
                            u[1][2]=awy*awz*(1.0-acos)-awx*asin;
                            u[2][0]=awx*awz*(1.0-acos)-awy*asin;
                            u[2][1]=awy*awz*(1.0-acos)+awx*asin;
                            u[2][2]=acos+awz*awz*(1.0-acos);            
                            // Rotation points
                    //      int rand_point = rand()%(pnumx);
                    //      int rand_point = rand()%(pnum);
                    //      axyz[0]=pointsx[rand_point].ca[0];
                    //      axyz[1]=pointsx[rand_point].ca[1];
                    //      axyz[2]=pointsx[rand_point].ca[2];
                    //        vector<float> coor_pdb(3,0.0);   
                    //        for(int j=rand_a;j<=rand_b;j++)
                    //        {
                    //            coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
                    //            coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
                    //            coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];             
                    //        }
                    //        coor_pdb[0] = coor_pdb[0]/abs(rand_b-rand_a+1);
                    //        coor_pdb[1] = coor_pdb[1]/abs(rand_b-rand_a+1);
                    //        coor_pdb[2] = coor_pdb[2]/abs(rand_b-rand_a+1);                                                   
                            int rand_point = rand_a;
                    //        axyz[0]= coor_pdb[0];  
                    //        axyz[1]= coor_pdb[1];
                    //        axyz[2]= coor_pdb[2];               
                            axyz[0] =pointsBx[3*rand_point+1].x_[0]; // CA atom
                            axyz[1] =pointsBx[3*rand_point+1].x_[1];
                            axyz[2] =pointsBx[3*rand_point+1].x_[2];
                            int tmp3=0;
                            vector<poseCoord > fin_mat;
                        //    fin_mat = tmp_mat;
                            vector<poseCoord>().swap(fin_mat);
                            fin_mat = fin_maty;
                      //      for(int j=rand_a+2;j<=rand_b-2;j++)
                        //    {
                        //        fin_mat.push_back(pointsBx[3*j+0]);
                        //        fin_mat.push_back(pointsBx[3*j+1]);
                        //        fin_mat.push_back(pointsBx[3*j+2]);               
                        //        tmp3=tmp3+1;
                        //    }        
                        //    old_dE = old_clash + olddis_p_m;                      
                            vector<poseCoord > fin_matx;
                            vector<poseCoord >().swap(fin_matx);
                            fin_matx=fin_mat;
                            tmp3=0;
                            float rand_mov0=(randf0and1()*2.0-1.0);
                            if(rand_mov0>0)
                            {
                                rand_point = rand_a;
                                axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
                                axyz[1]=pointsBx[3*rand_point+1].x_[1];
                                axyz[2]=pointsBx[3*rand_point+1].x_[2];
                                for(int j=3*rand_a+1;j<=3*rand_b+2;j++) // not last N
                                {
                                    fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                                    fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                                    fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                                    tmp3=tmp3+1;
                                }                                                        
                            } else 
                            {
                                rand_point = rand_b;
                                axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
                                axyz[1]=pointsBx[3*rand_point+1].x_[1];
                                axyz[2]=pointsBx[3*rand_point+1].x_[2];
                                for(int j=3*rand_a;j<=3*rand_b+1;j++) // not last N
                                {
                                    fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                                    fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                                    fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                                    tmp3=tmp3+1;
                                }                            
                            }
                            vector<poseCoord >().swap(fin_matx);
                            int numseq=pnum;
                        //    if((rand_b-rand_a-3)<2) continue;
                        //    if(rand_a==0) continue;
                            vector<point3f> decstr(numseq);
                            for(int j=0;j<pnum;j++)
                            {
                                decstr[j].ss2=numtoss(bb[j].sst);
                                decstr[j].ssm=numtoss(bb[j].sst);
                                decstr[j].stype=numtoss(bb[j].sst);
                                decstr[j].x=fin_mat[3*j+1].x_[0];
                                decstr[j].y=fin_mat[3*j+1].x_[1];
                                decstr[j].z=fin_mat[3*j+1].x_[2];
                                decstr[j].ptn.x=fin_mat[3*j+0].x_[0];
                                decstr[j].ptn.y=fin_mat[3*j+0].x_[1];
                                decstr[j].ptn.z=fin_mat[3*j+0].x_[2];
                                decstr[j].ptc.x=fin_mat[3*j+2].x_[0];
                                decstr[j].ptc.y=fin_mat[3*j+2].x_[1];
                                decstr[j].ptc.z=fin_mat[3*j+2].x_[2];
                //                tmp_tt=tmp_tt+1;
                            }
                            if(rand_mov0>0)
                            {
                                if((rand_b2)>(pnum-1)||(rand_b-1)<1) continue;
                                mcfragsweepLMP2(decstr,numseq,rand_b-1,rand_b2);
                            } else
                            {
                                if((rand_a2)<1||(rand_a+2)>(pnum-1)) continue;
                                mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+2);
                            }
                            
                            tmp_tt=0;
                            int tmp_ttx=0;
                            for(int j=0;j<pnum;j++)
                            {
                                fin_maty[3*j+1].x_[0]=decstr[j].x;
                                fin_maty[3*j+1].x_[1]=decstr[j].y;
                                fin_maty[3*j+1].x_[2]=decstr[j].z;
                                fin_maty[3*j+0].x_[0]=decstr[j].ptn.x;
                                fin_maty[3*j+0].x_[1]=decstr[j].ptn.y;
                                fin_maty[3*j+0].x_[2]=decstr[j].ptn.z;
                                fin_maty[3*j+2].x_[0]=decstr[j].ptc.x;
                                fin_maty[3*j+2].x_[1]=decstr[j].ptc.y;
                                fin_maty[3*j+2].x_[2]=decstr[j].ptc.z;
                                tmp_tt=tmp_tt+1;
                            }     
                            vector<point3f>().swap(decstr);  
                            NN2TF=true;              
                        }     
                       if(rand_movement==4)
                        {
                    /*        bool flagok0=false;
                            for(int j=rand_a;j<=rand_b;j++)
                            {
                                if(numtoss(bb[j].sst)!='C') flagok0 = true;
                            }
                            if(flagok0) continue; */
                //            if((rand_b+1)>(pnum-1)) continue;
                            float rand_mov0=(randf0and1()*2.0-1.0);
                            if(rand_mov0>=0.0)
                            {
                                vector<float> sftv(3,0.0);
                                sftv[0] = fin_maty[3*rand_b+1].x_[0] -fin_maty[3*rand_a+1].x_[0];
                                sftv[1] = fin_maty[3*rand_b+1].x_[1] -fin_maty[3*rand_a+1].x_[1];
                                sftv[2] = fin_maty[3*rand_b+1].x_[2] -fin_maty[3*rand_a+1].x_[2];
                                sftv[0] = sftv[0]/float(3*rand_b-3*rand_a+1);
                                sftv[1] = sftv[1]/float(3*rand_b-3*rand_a+1);
                                sftv[2] = sftv[2]/float(3*rand_b-3*rand_a+1);
    //                            sftv[0] = (2.0*randf0and1()-1.0);
    //                            sftv[1] = (2.0*randf0and1()-1.0);
    //                            sftv[2] = (2.0*randf0and1()-1.0);
                                for(int j=rand_a;j<=rand_b;j++)
                                {
                                    fin_maty[3*j+0].x_[0] = fin_maty[3*(j)+0].x_[0] + sftv[0];
                                    fin_maty[3*j+0].x_[1] = fin_maty[3*(j)+0].x_[1] + sftv[1];
                                    fin_maty[3*j+0].x_[2] = fin_maty[3*(j)+0].x_[2] + sftv[2];
                                    fin_maty[3*j+1].x_[0] = fin_maty[3*(j)+1].x_[0] + sftv[0];
                                    fin_maty[3*j+1].x_[1] = fin_maty[3*(j)+1].x_[1] + sftv[1];
                                    fin_maty[3*j+1].x_[2] = fin_maty[3*(j)+1].x_[2] + sftv[2];
                                    fin_maty[3*j+2].x_[0] = fin_maty[3*(j)+2].x_[0] + sftv[0];
                                    fin_maty[3*j+2].x_[1] = fin_maty[3*(j)+2].x_[1] + sftv[1];
                                    fin_maty[3*j+2].x_[2] = fin_maty[3*(j)+2].x_[2] + sftv[2];
                                }
                                int numseq=pnum;
                            //    if((rand_b-rand_a-3)<2) continue;
                            //    if(rand_a==0) continue;
                                vector<point3f> decstr(numseq);

                                for(int j=0;j<pnum;j++)
                                {
                                    decstr[j].ss2=numtoss(bb[j].sst);
                                    decstr[j].ssm=numtoss(bb[j].sst);
                                    decstr[j].stype=numtoss(bb[j].sst);
                                    decstr[j].x=fin_maty[3*j+1].x_[0];
                                    decstr[j].y=fin_maty[3*j+1].x_[1];
                                    decstr[j].z=fin_maty[3*j+1].x_[2];
                                    decstr[j].ptn.x=fin_maty[3*j+0].x_[0];
                                    decstr[j].ptn.y=fin_maty[3*j+0].x_[1];
                                    decstr[j].ptn.z=fin_maty[3*j+0].x_[2];
                                    decstr[j].ptc.x=fin_maty[3*j+2].x_[0];
                                    decstr[j].ptc.y=fin_maty[3*j+2].x_[1];
                                    decstr[j].ptc.z=fin_maty[3*j+2].x_[2];
                    //                tmp_tt=tmp_tt+1;
                                }
                    ///            int rand_a_a = rand_a-2;
                    //            int rand_a_b = rand_a+2;
                    //            if(rand_a+2>pnum-1) rand_a_b = pnum-1;
                    //            if((rand_a-2)<1) rand_a_a = 1;
                    //            mcfragsweepLMP2(decstr,numseq,rand_a_a,rand_a_b);
                                if(rand_a2<1) continue;
                                if(rand_a+2>(pnum-1)) continue;
                                if(rand_b-2<1) continue;
                                if(rand_b2>pnum-1) continue;                                
                                mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+2);
                    //            int rand_b_a = rand_b-2;
                    //            int rand_b_b = rand_b+2;
                    //            if((rand_b-2)<1) rand_b_a=1;
                    //            if(rand_b+2>pnum-1) rand_b_b = pnum-1;
                    //            mcfragsweepLMP2(decstr,numseq,rand_b_a,rand_b_b);
                                mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b2);
                                for(int j=0;j<pnum;j++)
                                {
                                    fin_maty[3*j+1].x_[0]=decstr[j].x;
                                    fin_maty[3*j+1].x_[1]=decstr[j].y;
                                    fin_maty[3*j+1].x_[2]=decstr[j].z;
                                    fin_maty[3*j+0].x_[0]=decstr[j].ptn.x;
                                    fin_maty[3*j+0].x_[1]=decstr[j].ptn.y;
                                    fin_maty[3*j+0].x_[2]=decstr[j].ptn.z;
                                    fin_maty[3*j+2].x_[0]=decstr[j].ptc.x;
                                    fin_maty[3*j+2].x_[1]=decstr[j].ptc.y;
                                    fin_maty[3*j+2].x_[2]=decstr[j].ptc.z;
                                }     
                                vector<point3f>().swap(decstr);
                            } else
                            {
                                vector<float> sftv(3,0.0);
                                sftv[0] = -fin_maty[3*rand_b+1].x_[0] + fin_maty[3*rand_a+1].x_[0];
                                sftv[1] = -fin_maty[3*rand_b+1].x_[1] + fin_maty[3*rand_a+1].x_[1];
                                sftv[2] = -fin_maty[3*rand_b+1].x_[2] + fin_maty[3*rand_a+1].x_[2];
                                sftv[0] = sftv[0]/float(3*rand_b-3*rand_a+1);
                                sftv[1] = sftv[1]/float(3*rand_b-3*rand_a+1);
                                sftv[2] = sftv[2]/float(3*rand_b-3*rand_a+1);
    //                            sftv[0] = (2.0*randf0and1()-1.0);
    //                            sftv[1] = (2.0*randf0and1()-1.0);
    //                            sftv[2] = (2.0*randf0and1()-1.0);
                                for(int j=rand_a;j<=rand_b;j++)
                                {
                                    fin_maty[3*j+0].x_[0] = fin_maty[3*(j)+0].x_[0] + sftv[0];
                                    fin_maty[3*j+0].x_[1] = fin_maty[3*(j)+0].x_[1] + sftv[1];
                                    fin_maty[3*j+0].x_[2] = fin_maty[3*(j)+0].x_[2] + sftv[2];
                                    fin_maty[3*j+1].x_[0] = fin_maty[3*(j)+1].x_[0] + sftv[0];
                                    fin_maty[3*j+1].x_[1] = fin_maty[3*(j)+1].x_[1] + sftv[1];
                                    fin_maty[3*j+1].x_[2] = fin_maty[3*(j)+1].x_[2] + sftv[2];
                                    fin_maty[3*j+2].x_[0] = fin_maty[3*(j)+2].x_[0] + sftv[0];
                                    fin_maty[3*j+2].x_[1] = fin_maty[3*(j)+2].x_[1] + sftv[1];
                                    fin_maty[3*j+2].x_[2] = fin_maty[3*(j)+2].x_[2] + sftv[2];
                                }
                                int numseq=pnum;
                            //    if((rand_b-rand_a-3)<2) continue;
                            //    if(rand_a==0) continue;
                                vector<point3f> decstr(numseq);

                                for(int j=0;j<pnum;j++)
                                {
                                    decstr[j].ss2=numtoss(bb[j].sst);
                                    decstr[j].ssm=numtoss(bb[j].sst);
                                    decstr[j].stype=numtoss(bb[j].sst);
                                    decstr[j].x=fin_maty[3*j+1].x_[0];
                                    decstr[j].y=fin_maty[3*j+1].x_[1];
                                    decstr[j].z=fin_maty[3*j+1].x_[2];
                                    decstr[j].ptn.x=fin_maty[3*j+0].x_[0];
                                    decstr[j].ptn.y=fin_maty[3*j+0].x_[1];
                                    decstr[j].ptn.z=fin_maty[3*j+0].x_[2];
                                    decstr[j].ptc.x=fin_maty[3*j+2].x_[0];
                                    decstr[j].ptc.y=fin_maty[3*j+2].x_[1];
                                    decstr[j].ptc.z=fin_maty[3*j+2].x_[2];
                    //                tmp_tt=tmp_tt+1;
                                }
                        //        int rand_a_a = rand_a-2;
                        //        int rand_a_b = rand_a+2;
                        //        if(rand_a+2>pnum-1) rand_a_b = pnum-1;
                        //        if((rand_a-2)<1) rand_a_a = 1;
                        //        mcfragsweepLMP2(decstr,numseq,rand_a_a,rand_a_b);
                                if(rand_a2<1) continue;
                                if(rand_a+2>(pnum-1)) continue;
                                if(rand_b-2<1) continue;
                                if(rand_b2>pnum-1) continue;                                
                                mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+2);
                        //        int rand_b_a = rand_b-2;
                        //        int rand_b_b = rand_b+2;
                        //        if((rand_b-2)<1) rand_b_a=1;
                        //        if(rand_b+2>pnum-1) rand_b_b = pnum-1;
                        //        mcfragsweepLMP2(decstr,numseq,rand_b_a,rand_b_b);
                                mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b2);
                                for(int j=0;j<pnum;j++)
                                {
                                    fin_maty[3*j+1].x_[0]=decstr[j].x;
                                    fin_maty[3*j+1].x_[1]=decstr[j].y;
                                    fin_maty[3*j+1].x_[2]=decstr[j].z;
                                    fin_maty[3*j+0].x_[0]=decstr[j].ptn.x;
                                    fin_maty[3*j+0].x_[1]=decstr[j].ptn.y;
                                    fin_maty[3*j+0].x_[2]=decstr[j].ptn.z;
                                    fin_maty[3*j+2].x_[0]=decstr[j].ptc.x;
                                    fin_maty[3*j+2].x_[1]=decstr[j].ptc.y;
                                    fin_maty[3*j+2].x_[2]=decstr[j].ptc.z;
                                }     
                                vector<point3f>().swap(decstr);
                            }
                            NN4TF=true;
                        } 
                    if(rand_movement==5)
                    {
                /*        bool flagok0=false;
                        for(int j=rand_a;j<=rand_b;j++)
                        {
                            if(numtoss(bb[j].sst)!='C') flagok0 = true;
                        }
                        if(flagok0) continue; */
            //            if((rand_b+1)>(pnum-1)) continue;
                        vector<float> sftv(3,0.0);
                        sftv[0] = (2.0*randf0and1()-1.0);
                        sftv[1] = (2.0*randf0and1()-1.0);
                        sftv[2] = (2.0*randf0and1()-1.0);
                        for(int j=rand_a;j<=rand_b;j++)
                        {
                            fin_maty[3*j+0].x_[0] = fin_maty[3*(j)+0].x_[0] + sftv[0];
                            fin_maty[3*j+0].x_[1] = fin_maty[3*(j)+0].x_[1] + sftv[1];
                            fin_maty[3*j+0].x_[2] = fin_maty[3*(j)+0].x_[2] + sftv[2];
                            fin_maty[3*j+1].x_[0] = fin_maty[3*(j)+1].x_[0] + sftv[0];
                            fin_maty[3*j+1].x_[1] = fin_maty[3*(j)+1].x_[1] + sftv[1];
                            fin_maty[3*j+1].x_[2] = fin_maty[3*(j)+1].x_[2] + sftv[2];
                            fin_maty[3*j+2].x_[0] = fin_maty[3*(j)+2].x_[0] + sftv[0];
                            fin_maty[3*j+2].x_[1] = fin_maty[3*(j)+2].x_[1] + sftv[1];
                            fin_maty[3*j+2].x_[2] = fin_maty[3*(j)+2].x_[2] + sftv[2];
                        }
                        int numseq=pnum;
                    //    if((rand_b-rand_a-3)<2) continue;
                    //    if(rand_a==0) continue;
                        vector<point3f> decstr(numseq);

                        for(int j=0;j<pnum;j++)
                        {
                            decstr[j].ss2=numtoss(bb[j].sst);
                            decstr[j].ssm=numtoss(bb[j].sst);
                            decstr[j].stype=numtoss(bb[j].sst);
                            decstr[j].x=fin_maty[3*j+1].x_[0];
                            decstr[j].y=fin_maty[3*j+1].x_[1];
                            decstr[j].z=fin_maty[3*j+1].x_[2];
                            decstr[j].ptn.x=fin_maty[3*j+0].x_[0];
                            decstr[j].ptn.y=fin_maty[3*j+0].x_[1];
                            decstr[j].ptn.z=fin_maty[3*j+0].x_[2];
                            decstr[j].ptc.x=fin_maty[3*j+2].x_[0];
                            decstr[j].ptc.y=fin_maty[3*j+2].x_[1];
                            decstr[j].ptc.z=fin_maty[3*j+2].x_[2];
            //                tmp_tt=tmp_tt+1;
                        }
                //        if(randtx>(pnum-7))
                //        {
                //            if((rand_a-2)<1||(rand_a+2)>(pnum-1)) continue;
                            if(rand_a2<1) continue;
                            if(rand_a+2>(pnum-1)) continue;
                            if(rand_b-2<1) continue;
                            if(rand_b2>pnum-1) continue;                        
                            mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+2);
                //        } else
                //        {
                //            if((rand_b+2)>(pnum-1)||(rand_b-2)<1) continue;
                            mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b2); 
                //        }
                        for(int j=0;j<pnum;j++)
                        {
                            fin_maty[3*j+1].x_[0]=decstr[j].x;
                            fin_maty[3*j+1].x_[1]=decstr[j].y;
                            fin_maty[3*j+1].x_[2]=decstr[j].z;
                            fin_maty[3*j+0].x_[0]=decstr[j].ptn.x;
                            fin_maty[3*j+0].x_[1]=decstr[j].ptn.y;
                            fin_maty[3*j+0].x_[2]=decstr[j].ptn.z;
                            fin_maty[3*j+2].x_[0]=decstr[j].ptc.x;
                            fin_maty[3*j+2].x_[1]=decstr[j].ptc.y;
                            fin_maty[3*j+2].x_[2]=decstr[j].ptc.z;
                        }     
                        vector<point3f>().swap(decstr);  
                        NN3TF = true;              
                    }                    
                }

                    dis_p=0.0;
                 //   tp_clashx= false;
                    new_vwd = 0.0;
                    new_clash = 0.0;
                    for(int jj=3*rand_a2;jj<=3*rand_b2+2;jj++)
                    {
                        if(fin_maty[jj].elt_=="CA")
                        {
                //            float VR1= GetVdwRadius(tmp_mat[jj].elt_);
                            string VR1 = fin_maty[jj].elt_;
                            for(int t=0;t<fin_maty.size();t++)
                            {
                                if(fin_maty[t].elt_=="CA")
                                {
                    //                if((jj+rand_a)==t) continue;
                            //        if(t>=(3*rand_a)&& t<=(3*rand_b+2)) continue;
                                    if(abs(t-jj)<2) continue;
                    //                if(t>3*rand_a && t<(3*rand_a+tmp_mat.size())) continue;
                    //                float VR2= GetVdwRadius(pointsBx[t].elt_);
                                    string VR2 = fin_maty[t].elt_ ;
                                    dis_p = Distance_point(fin_maty[jj].x_,fin_maty[t].x_);
                                    //if(dis_p < (VR1+VR2))
                            //        if(dis_p < 3.70)
                            //        {
                            //            new_clash = new_clash + 1.0/ sqrt(dis_p + 0.0001) ;
                     //                   old_clash = old_clash + exp( 3.75-dis_p) ;
                            //        }
                                      new_clash = new_clash + GetVdwEgCG(VR1,VR2,dis_p);                        
                    //                new_vwd = new_vwd + GetVdwEg(VR1,VR2,dis_p);
                    //                new_vwd = new_vwd + GetVdwEgC(VR1,VR2,dis_p);
                    //                if((dis_p)<2.0) tp_clashx = true;
                    //                if(tp_clashx)
                    //                {
                    //                    cout<<"jj,t: "<<jj+rand_a<<" "<<t<<" dis: "<<VR1+VR2<<" "<<dis_p<<endl;
                    //                    tp_clashx=false;
                    //                }
                //                    cout<< "VR1+VR2,dis_p: "<<VR1+VR2<<" "<<dis_p<<endl;
                                }
                            }
            //                if(tp_clashx) break;
                        }
                    }
            //            cout<<"tp_clashx: " <<tp_clashx<<endl;
            //            if(tp_clashx) continue;  
            //            cout<<" 3 "<<endl;                 
                    ObjexxFCL::FArray3D< float > rhoC2;
                    ObjexxFCL::FArray3D< float > inv_rho_mask2;
                    rhoC2.dimension(pDenx , pDeny , pDenz);
                    inv_rho_mask2.dimension(pDenx , pDeny , pDenz);
                        for ( int t=0; t<pDenxyz; ++t ) {
                            rhoC2[t]=0.0;
                            inv_rho_mask2[t]=1.0;
                        }
                    vector<float> del_ij(3,0.0);
                    vector<float> atm_j(3,0.0);             
                    del_ij=vector<float>(3,0.0);
                    atm_j=vector<float>(3,0.0);
                    int tm_i1 =0;
                    vector<float> atm_idy;
                    atm_idy=vector<float>(3,0.0);
                    bond0max=-1000.0;bond1max=-1000.0;bond2max=-1000.0;
                    bond0min=1000.0;bond1min=1000.0;bond2min=1000.0;

                    for(int j=3*rand_a2;j<=3*rand_b2+2;j++)
                    {
                        if(fin_maty[j].elt_ == "CA")
                        {
                            vector<float> cartX1;
                            vector<float> fracX1; 
                        //    elt_i = tmp_mat[3*tm_i1+1].elt_;
                            elt_i = fin_maty[j].elt_;
                            elt_i = elt_i[0];
                            OneGaussianScattering sig_j = get_A( elt_i );
                            k = sig_j.k( theDensityMapx.effectiveB );
                            C = sig_j.C( k );
                            if ( C < 1e-6 ) continue;  

                    //        cartX1 = tmp_mat[3*tm_i1+1].x_;
                            cartX1 = fin_maty[j].x_;
                    //        tm_i1 = tm_i1 + 1;
                            MatrixTimesTransVector(theDensityMapx.c2f,cartX1,fracX1);
                            atm_idy[0] = pos_mod (double(fracX1[0]*pgrid[0] - porigin[0] + 1) , (double)pgrid[0]);
                            atm_idy[1] = pos_mod (double(fracX1[1]*pgrid[1] - porigin[1] + 1) , (double)pgrid[1]);
                            atm_idy[2] = pos_mod (double(fracX1[2]*pgrid[2] - porigin[2] + 1) , (double)pgrid[2]);                 
                            for(int z=1;z<=pDenz;z++)
                            {
                                atm_j[2] = z;
                                del_ij[2] =(atm_idy[2]-atm_j[2])/pgrid[2];
                                if ( del_ij[2] > 0.5 ) del_ij[2]-=1.0;
                                if ( del_ij[2] < -0.5 ) del_ij[2]+=1.0;                    
                                del_ij[0] = del_ij[1] = 0.0;
                                vector<float> frac_tmpz;
                                MatrixTimesTransVector(theDensityMapx.f2c,del_ij,frac_tmpz);
                                if(square_len(frac_tmpz)> (ppadding+patom_m)*(ppadding+patom_m) ) continue; 
                                if(z>bond2max) bond2max = z; 
                                if(z<bond2min) bond2min = z;                
                                for(int y=1;y<=pDeny;y++)
                                {
                                    atm_j[1] = y;
                                    del_ij[1] = (atm_idy[1] - atm_j[1])/pgrid[1];
                                    // wrap-around??
                                    if ( del_ij[1] > 0.5 ) del_ij[1]-=1.0;
                                    if ( del_ij[1] < -0.5 ) del_ij[1]+=1.0;                        
                                    del_ij[0] = 0.0;
                                    vector<float> frac_tmpy;
                                    MatrixTimesTransVector(theDensityMapx.f2c,del_ij,frac_tmpy);
                                    if(square_len(frac_tmpy)> (ppadding+patom_m)*(ppadding+patom_m) ) continue;                                
                                    if(y<bond1min) bond1min = y;
                                    if(y>bond1max) bond1max = y;                  
                                    for(int x=1;x<=pDenx;x++)
                                    {
                                        atm_j[0] = x;
                                        del_ij[0] = (atm_idy[0] - atm_j[0])/pgrid[0];
                                        // wrap-around??
                                        if ( del_ij[0] > 0.5 ) del_ij[0]-=1.0;
                                        if ( del_ij[0] < -0.5 ) del_ij[0]+=1.0;                            
                                        vector<float> cart_del_ij2;
                                        MatrixTimesTransVector(theDensityMapx.f2c,del_ij,cart_del_ij2);
                                        float d2 = square_len(cart_del_ij2);
                                        if(d2 > (ppadding+patom_m)*(ppadding+patom_m) ) continue;
                                    
                                        float atm = C*exp(-k*d2);
                                        float sigmoid_msk = exp( d2 - (patom_m)*(patom_m)  );
                                    //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK));
                                    //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK));
                                        float inv_msk = 1.0/(1.0+sigmoid_msk);
                        //                rhoC2(x,y,z) += atm;
                        //                rhoC2(x,y,z) += atm;
                                        rhoC2(x,y,z) += atm;
                                        inv_rho_mask2(x,y,z) *= (1.0 - inv_msk);


                                        if(x>bond0max) bond0max = x;
                                        
                                        
                                        if(x<bond0min) bond0min = x;
                                        
                                        
            //                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
            //                                mask2(x,y,z) = 1.0; // problem?
            //                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
            //                                    float atm = C*exp(-k*d2);
            //                                    rhoC(x,y,z) += atm;
            //                                }
            //                            }
                                    }
                                }
                            }
                            tm_i1 = tm_i1 + 1;
                        }
                    }
            //        rhoC2 = rhoC0;
        //            max_rhc = 0.0;
        //            for ( int x=0; x<pDenxyz; ++x ) {
        //                if(rhoC2[x]>max_rhc) max_rhc = rhoC2[x] ;
        //            }     
            /*        for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
                //        cout<<theDensityMap.density[x]<<" ";
                        if(rhoC2[x]<(3.0/5.0)*max_rhc) rhoC2[x]=0; 
                    } */
        //            rhoC2 = rhoC0;                           
                    sumC_i2=0.0, sumO_i2=0.0, sumCO_i2=0.0, vol_i2=0.0, CC_i2=0.0;
                    sumO2_i2=0.0, sumC2_i2=0.0, varC_i2=0.0, varO_i2=0.0;
                    clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
            /*            for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
                //            tmp_den=0.0;
                //            tmp_rhc =0.0;
                //            if(theDensityMap.density[x] > (2.0/3.0)*max_den) tmp_den= theDensityMap.density[x];
                //            if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
                //            if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                            // fetch this point
                            clc_x2 = rhoC2[x];
                //            clc_x2 = tmp_rhc;
                            obs_x2 = theDensityMap.density[x];
                //            obs_x2 = tmp_den;
                            eps_x2 = 1-inv_rho_mask2[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                //            eps_x2 = 1.0;

                            // SMOOTHED
                            sumCO_i2 += eps_x2*clc_x2*obs_x2;
                            sumO_i2  += eps_x2*obs_x2;
                            sumO2_i2 += eps_x2*obs_x2*obs_x2;
                            sumC_i2  += eps_x2*clc_x2;
                            sumC2_i2 += eps_x2*clc_x2*clc_x2;
                            vol_i2   += eps_x2;
                        }
                        varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
                        varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
                        if ( varC_i2 == 0 || varO_i2 == 0 ) {
                            CC_i2 = 0;
                        } else {
                            CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
                        }   */
                    for(int x=bond0min;x<=bond0max;x++)
                    {
                       for(int y=bond1min;y<=bond1max;y++)
                        {
                            for(int z=bond2min;z<=bond2max;z++)
                            {
                   //            tmp_den=0.0;
                    //            tmp_rhc =0.0;
                    //            if(theDensityMap.density[x] > (2.0/3.0)*max_den) tmp_den= theDensityMap.density[x];
                    //            if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
                    //            if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                                // fetch this point
                                clc_x2 = rhoC2(x,y,z);
                    //            clc_x2 = tmp_rhc;
                                obs_x2 = theDensityMapx.density(x,y,z);
                    //            obs_x2 = tmp_den;
                                eps_x2 = 1-inv_rho_mask2(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                    //            eps_x2 = 1.0;

                                // SMOOTHED
                                sumCO_i2 += eps_x2*clc_x2*obs_x2;
                                sumO_i2  += eps_x2*obs_x2;
                                sumO2_i2 += eps_x2*obs_x2*obs_x2;
                                sumC_i2  += eps_x2*clc_x2;
                                sumC2_i2 += eps_x2*clc_x2*clc_x2;
                                vol_i2   += eps_x2;
                            }
                        }
                    }          
                    varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
                    varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 );
                    if ( varC_i2 == 0 || varO_i2 == 0 ) {
                        CC_i2 = 0;
                    } else {
                        CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
                    }

                        new_CC1 =1.0-CC_i2; 

            /*        vector<point3f> decstrbondangx(pnum);
                    for(int j=0;j<pnum;j++)
                    {
                        decstrbondangx[j].ss2=numtoss(bb[j].sst);
                        decstrbondangx[j].x=fin_maty[3*j+1].x_[0];
                        decstrbondangx[j].y=fin_maty[3*j+1].x_[1];
                        decstrbondangx[j].z=fin_maty[3*j+1].x_[2];
                        decstrbondangx[j].ptn.x=fin_maty[3*j+0].x_[0];
                        decstrbondangx[j].ptn.y=fin_maty[3*j+0].x_[1];
                        decstrbondangx[j].ptn.z=fin_maty[3*j+0].x_[2];
                        decstrbondangx[j].ptc.x=fin_maty[3*j+2].x_[0];
                        decstrbondangx[j].ptc.y=fin_maty[3*j+2].x_[1];
                        decstrbondangx[j].ptc.z=fin_maty[3*j+2].x_[2];
                        string resng = pointsrenm[j];
                        char resny = resng[0];
                        decstrbondangx[j].iaa = aminoid(resny);
            //                decstr[i].iaa=pp.aminoid(decstr[i].aaa);                
            //                tmp_tt=tmp_tt+1;
                    }     
                    for(int j=0;j<pnum;j++)
                    {
                        if(j==0)
                        {
                            decstrbondangx[j].ang[0]=0.0;
                            decstrbondangx[j].ang[1]=0.0;
                            decstrbondangx[j].ang[2]=0.0;
                        }               
                        point3d p12 = setv(decstrbondangx[j-1].x-decstrbondangx[j-1].ptc.x,decstrbondangx[j-1].y-decstrbondangx[j-1].ptc.y,decstrbondangx[j-1].z-decstrbondangx[j-1].ptc.z);
                        point3d p23 = setv(decstrbondangx[j].ptn.x-decstrbondangx[j-1].ptc.x,decstrbondangx[j].ptn.y-decstrbondangx[j-1].ptc.y,decstrbondangx[j].ptn.z-decstrbondangx[j-1].ptc.z);
                        decstrbondangx[j].ang[0]=angv(p12,p23)*degrad;
                        p12 = setv(decstrbondangx[j-1].ptc.x-decstrbondangx[j].ptn.x,decstrbondangx[j-1].ptc.y-decstrbondangx[j].ptn.y,decstrbondangx[j-1].ptc.z-decstrbondangx[j].ptn.z);
                        p23 = setv(decstrbondangx[j].x-decstrbondangx[j].ptn.x,decstrbondangx[j].y-decstrbondangx[j].ptn.y,decstrbondangx[j].z-decstrbondangx[j].ptn.z);
                        decstrbondangx[j].ang[1]=angv(p12,p23)*degrad;
                        p12 = setv(decstrbondangx[j].ptn.x-decstrbondangx[j].x,decstrbondangx[j].ptn.y-decstrbondangx[j].y,decstrbondangx[j].ptn.z-decstrbondangx[j].z);
                        p23 = setv(decstrbondangx[j].ptc.x-decstrbondangx[j].x,decstrbondangx[j].ptc.y-decstrbondangx[j].y,decstrbondangx[j].ptc.z-decstrbondangx[j].z);
                        decstrbondangx[j].ang[2]= angv(p12,p23)*degrad;
            //                tmp_tt=tmp_tt+1;
                    }
                    float new_bondangenergy =0.0;
                    Eenergybondangle(decstrbondangx,pnum,new_bondangenergy);
                    cout<< "new_bondangenergy: "<< new_bondangenergy/pnum<<endl;  */
                        float new_Ehbond=0.0;
                //        point3f *new_decstr;
                     //   if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E')
                        {                       
                    /*        tmp_tt=0;
                            for(int j=0;j<pnum;j++)
                            {
                                bby[j].indn=3*j+0;
                                bby[j].indca=3*j+1;
                                bby[j].indc=3*j+2; 
                                string resng = pointsrenm[j];
                                char resny = resng[0];
                                bby[j].resid = aminoid(resny);                                    
                                tmp_tt=tmp_tt+1;
                            }     
                        //    calcsse2(bby, pnum, tmp_mat); 
                            calcssennhoc(bby,pnum,tmp_mat);   */
                            vector<boneinfo> bbg(pnum);
                            for(int j=0;j<pnum;j++)
                            {
                                bbg[j].indn = 3*j+0;
                                bbg[j].indca = 3*j+1;
                                bbg[j].indc = 3*j+2; 
                                string resng = pointsrenm[j];
                                char resny = resng[0];
                                bbg[j].resid = aminoid(resny);                                    
                            }     
                        //    calcsse2(bbg, pnum, tmp_mat); 
                            calcssennhoc(bbg,pnum,fin_maty);                

        //                if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E')
        //                {
                            int numseq=pnum;
                            vector<point3f> new_decstr(numseq);             
                            tmp_tt=0;
                            for(int j=0;j<pnum;j++)
                            {                     
                            //    new_decstr[tmp_tt].ss2=numtoss(bby[j].sst);
                                new_decstr[j].ss2=numtoss(bbg[j].sst);
                                new_decstr[j].ssm=numtoss(bbg[j].sst);
                                new_decstr[j].x=fin_maty[3*tmp_tt+1].x_[0];
                                new_decstr[j].y=fin_maty[3*tmp_tt+1].x_[1];
                                new_decstr[j].z=fin_maty[3*tmp_tt+1].x_[2];
                                new_decstr[j].ptn.x=fin_maty[3*tmp_tt+0].x_[0];
                                new_decstr[j].ptn.y=fin_maty[3*tmp_tt+0].x_[1];
                                new_decstr[j].ptn.z=fin_maty[3*tmp_tt+0].x_[2];
                                new_decstr[j].ptc.x=fin_maty[3*tmp_tt+2].x_[0];
                                new_decstr[j].ptc.y=fin_maty[3*tmp_tt+2].x_[1];
                                new_decstr[j].ptc.z=fin_maty[3*tmp_tt+2].x_[2];
                                tmp_tt=tmp_tt+1;
                            }                  
                        //    new_Ehbond=energyhbondcanc(new_decstr,numseq);
                            new_Ehbond = energyhbondnhoc2(new_decstr,numseq);                           
                            vector<point3f>().swap(new_decstr);
                        }                    
            /*        if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E')
                    {                       
                        tmp_tt=0;
                        for(int j=0;j<pnum;j++)
                        {
                            bby[j].indn=3*j+0;
                            bby[j].indca=3*j+1;
                            bby[j].indc=3*j+2;
                            string resng = pointsrenm[j];
                            char resny = resng[0];
                            bb[j].resid = aminoid(resny);                                 
                            tmp_tt=tmp_tt+1;
                        }     
                    //    calcsse2(bby, pnum, tmp_mat); 
                        calcssennhoc(bby,pnum,tmp_mat);   
                        vector<boneinfo> bbg(pnum);
                        for(int j=0;j<pnum;j++)
                        {
                            bbg[j].indn = 3*j+0;
                            bbg[j].indca = 3*j+1;
                            bbg[j].indc = 3*j+2; 
                            string resng = pointsrenm[j];
                            char resny = resng[0];
                            bbg[j].resid = aminoid(resny);                               
                        }     
                    //    calcsse2(bbg, pnum, fin_maty); 
                        calcssennhoc(bbg,pnum,fin_maty);                    

    //                if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E')
    //                {
                        int numseq=pnum;
                        vector<point3f> new_decstr(numseq);             
                        tmp_tt=0;
                        for(int j=0;j<pnum;j++)
                        {                     
                        //    new_decstr[tmp_tt].ss2=numtoss(bby[j].sst);
                            new_decstr[j].ss2=numtoss(bbg[j].sst);
                            new_decstr[j].x=fin_maty[3*tmp_tt+1].x_[0];
                            new_decstr[j].y=fin_maty[3*tmp_tt+1].x_[1];
                            new_decstr[j].z=fin_maty[3*tmp_tt+1].x_[2];
                            new_decstr[j].ptn.x=fin_maty[3*tmp_tt+0].x_[0];
                            new_decstr[j].ptn.y=fin_maty[3*tmp_tt+0].x_[1];
                            new_decstr[j].ptn.z=fin_maty[3*tmp_tt+0].x_[2];
                            new_decstr[j].ptc.x=fin_maty[3*tmp_tt+2].x_[0];
                            new_decstr[j].ptc.y=fin_maty[3*tmp_tt+2].x_[1];
                            new_decstr[j].ptc.z=fin_maty[3*tmp_tt+2].x_[2];
                            tmp_tt=tmp_tt+1;
                        }                  
                    //    new_Ehbond=energyhbondcanc(new_decstr,numseq);
                        new_Ehbond=energyhbondnhoc2(new_decstr,numseq);
                        vector<point3f>().swap(new_decstr);
                    }   */                 


                //        cout<< "new_CC: "<<new_CC<<" ";
                    //    new_dE = 30*new_CC + new_vwd/((pnumx-1)*(rand_b-rand_a+1));
                    //    new_dE = 0.8*new_CC + 0.5*new_clash/abs(rand_b-rand_a+1) + 0.5*new_Ehbond/abs(rand_b-rand_a+1);
                    //    new_dE = 200.0*new_CC1 + 10.0*new_clash/abs(rand_b2-rand_a2+1) + 200.0*(new_Ehbond/abs(pnum))/(1.0-(new_Ehbond/abs(pnum)));//+new_bondangenergy/pnum;
                        float Eres=0.0;
                        for(int j=3*rand_a2;j<=3*rand_b2+2;j++)
                        {
                            for(int js=j-3;js<=j+3;js++)
                            {
                                if(js==j) continue;
                                if(js<0) continue;
                                if(js>(pnum-1)) continue;
                                float dx0 = Distance_point(fin_maty[j].x_,fin_maty[js].x_);
                                float dx1 = Distance_point(pointsBx0[j].x_,pointsBx0[js].x_);
                                Eres = Eres + (dx0-dx1)*(dx0-dx1);
                            }
                        }                        
                        Eres = Eres/float(rand_b2-rand_a2+1);

                        vector<point3f>().swap(decstrbondlen);
                        decstrbondlen = vector<point3f>(pnum);
                        for(int j=0;j<pnum;j++)
                        {
                            decstrbondlen[j].ss2=numtoss(bb[j].sst);
                            decstrbondlen[j].ssm=numtoss(bb[j].sst);
                            decstrbondlen[j].stype=numtoss(bb[j].sst);
                            decstrbondlen[j].x=fin_maty[3*j+1].x_[0];
                            decstrbondlen[j].y=fin_maty[3*j+1].x_[1];
                            decstrbondlen[j].z=fin_maty[3*j+1].x_[2];
                            decstrbondlen[j].ptn.x=fin_maty[3*j+0].x_[0];
                            decstrbondlen[j].ptn.y=fin_maty[3*j+0].x_[1];
                            decstrbondlen[j].ptn.z=fin_maty[3*j+0].x_[2];
                            decstrbondlen[j].ptc.x=fin_maty[3*j+2].x_[0];
                            decstrbondlen[j].ptc.y=fin_maty[3*j+2].x_[1];
                            decstrbondlen[j].ptc.z=fin_maty[3*j+2].x_[2];
                            string resng = pointsrenm[j];
                            char resny = resng[0];
                            decstrbondlen[j].iaa = aminoid(resny);
            //                decstr[i].iaa=pp.aminoid(decstr[i].aaa);                
            //                tmp_tt=tmp_tt+1;
                        } 
                        str2tor(decstrbondlen,pnum,3);

                        float E_distx = 0.0;
                        for(int iu=0;iu<N_dist;iu++)
                        {
                            int ik = dist_pre[iu].ires;
                            int jk = dist_pre[iu].jres;
                            float tmp_x = dist_pre[iu].dv;
                            float tmp_y = Distance_point(fin_maty[3*ik+1].x_,fin_maty[3*jk+1].x_);
                    //        float tmp_z=sqrt((tmp_x-tmp_y)*(tmp_x-tmp_y))-dist_pre[iu].ave;
                            float tmp_z=tmp_y-dist_pre[iu].ave;
                            E_distx = E_distx + (log(tmp_z*tmp_z+1))/sqrt(dist_pre[iu].stdx);

                        }    
                        E_distx=E_distx/float(pnum*(rand_b2-rand_a2+1));                        

                        float E_conx = 0.0;
                        for(int iu=0;iu<N_distg;iu++)
                        {
                            int ik = dist_map[iu].ires;
                            int jk = dist_map[iu].jres;
                            float tmp_x = dist_map[iu].dv;
                            float tmp_y = Distance_point(fin_maty[3*ik+1].x_,fin_maty[3*jk+1].x_);
                            E_conx = E_conx + sqrt((tmp_x-tmp_y)*(tmp_x-tmp_y));

                        }
                        E_conx=E_conx/float((rand_b2-rand_a2+1));
               /*         float new_fang=0.0;
                        float new_Ebondang=0.0;
                        new_Ebondang =(float) Eenergybondangle(decstrbondlen,pnum,new_fang);
                        new_Ebondang = 0.30*new_Ebondang + new_fang;      */                   
                    //    str2torp(decstrbondlen,pnum,0,pnum-1);  
                        float new_fene=0.0;
                        float new_Ebondlen=0.0;
                        new_Ebondlen = (float) energybondlength(decstrbondlen,pnum,new_fene);
                        new_Ebondlen = 0.5*new_Ebondlen + 20.0*new_fene;             

                        // torsion angle
                        float new_tor_fene=0.0;
                        float new_tor_E=0.0;
                        new_tor_E = energyrama(decstrbondlen,pnum,new_tor_fene,rlogduke,ramaduke); 
                        new_tor_E = 4.00*new_tor_E + new_tor_fene;          
                        vector<point3f>().swap(decstrbondlen);                          
                        new_dE = w1*new_CC1 + new_clash + new_Ebondlen + 10.0*Eres + w6*E_distx + new_tor_E + 1.0*E_conx + 0.5*new_Ehbond;// + new_Ebondang/1000.0;// + 200.0*(new_Ehbond/abs(pnum));
                
                //        cout<<"old_Ebondang, new_Ebondang: "<<old_Ebondang<<" "<<new_Ebondang<<endl; 
                    //    new_dE = new_CC1 + new_clash;
                    //    cout<<"new coor:  "<<tmp_mat[5].x_[0]<<" "<<tmp_mat[5].x_[1]<<" "<<tmp_mat[5].x_[2]<<endl;
                    //    cout<<"old: "<< 500.0*old_CC1<<" " <<old_clash<<endl;
                    //    cout<<"new: "<< 500.0*new_CC1<<" " <<new_clash<<endl;
                    //    cout<<"CC: "<<new_CC1<<" "<<old_CC1<<endl;
                    //   cout<<"dE_CC, dE_clash,new_Ebondlen,new_tor_E,new_Ehbond,E_distx,E_conx : "<< w1*(new_CC1-old_CC1)<<" "<<(new_clash-old_clash)<<" "<< new_Ebondlen-old_Ebondlen<<" "<<new_tor_E-old_tor_E <<" "<<0.5*(new_Ehbond-old_Ehbond)<<" "<<w6*(E_distx-E_dist)<<" "<<(E_conx -E_con) <<endl;
                //        cout<<"old new dE: "<<old_dE<<" "<<new_dE<<" "<<dE<<endl;
                        dE = new_dE - old_dE;
                //        cout<<"dE_CC, dE_clash,new_Ebondlen,new_tor_E,new_Ehbond,E_distx,E_conx : "<< w1*(new_CC1-old_CC1)<<" "<<(new_clash-old_clash)<<" "<< new_Ebondlen-old_Ebondlen<<" "<<new_tor_E-old_tor_E <<" "<<0.5*(new_Ehbond-old_Ehbond)<<" "<<w6*(E_distx-E_dist)<<" "<<(E_conx -E_con)<<" "<<dE <<endl;
                //        cout<<"dE: "<<dE<<endl;
                //        cout<<"dE_CC: "<<500.0*(new_CC1-old_CC1)<<endl;
                        if(new_dE <old_dE)
                        {
                //            cout<<"KT,new_dE: "<<KT<<" "<<new_dE<<endl;
                        //    cout<<"randx: "<<randtx<<endl;
                            int tmp_tmx=0;
                            for(int j=0;j<pnum;j++)
                            {
                                pointsBx[3*j+0] = fin_maty[3*j+0];
                                pointsBx[3*j+1] = fin_maty[3*j+1];
                                pointsBx[3*j+2] = fin_maty[3*j+2];
                            //    bb[j]=bby[tmp_tmx];
                                tmp_tmx = tmp_tmx +1;
                            }
                            if(new_dE<best_E)
                            {
                                best_E = new_dE;
                                vector<poseCoord>().swap(best_model);
                                best_model = fin_maty;
                            }  
                            for(int j=0;j<pnum;j++)
                            {
                                bb[j].indn = 3*j+0;
                                bb[j].indca = 3*j+1;
                                bb[j].indc = 3*j+2;  
                                string resng = pointsrenm[j];
                                char resny = resng[0];
                                bb[j].resid = aminoid(resny);                                       
                            }     
                        //    calcsse2(bb, pnum, tmp_mat); 
                            calcssennhoc(bb,pnum,fin_maty);                              
                                // calculate CC between the density map and ervery residue
                                vector<sssegment>().swap(segm);
                                if(calres)
                                {
                                    res_CC = vector<float>(pnum,2.0);
                                //    vector<float> res_CC3(pnum,0.0);
                                    for(int j=0;j<pnum;j++)
                                    {
                                        vector<poseCoord> sing_res(1);
                                        sing_res[0]= pointsBx[3*j+1];
                                        res_CC[j] = theDensityMap.matchposex(sing_res);
                                        res_CC3k[j]=res_CC[j];
                                //        res_CC3[j]=res_CC3k[j];

                                //        cout<<"j,CC: "<<j<<" "<<res_CC[j]<<endl; 
                                    }
                                    
                                //    nhg=0;
                                /*    for(int j=0;j<pnum;j++)
                                    {
                                        for(int jp=-1;jp<2;jp++)
                                        {
                                            int jjp = j+ jp;
                                            if(jjp<0) jjp=0;
                                            if(jjp>pnum-1) jjp=pnum-1;
                                            res_CC3[j] = res_CC3[j] + res_CC[jjp];
                                        }
                                        res_CC3k[j] = res_CC3[j]/3.0;
                        //                if(res_CC3k[j]<0.05) nhg=nhg+1;
                                    }  */
                        //            w1=w11*(float(pnum-nhg)/float(pnum))*(float(pnum-nhg)/float(pnum));
                        //            w6=w66*(float(pnum)/float(pnum-nhg))*(float(pnum)/float(pnum-nhg));               
                                    for(int j=0;j<pnum;j++)
                                    {
                                        int initx0=0;
                                        int initx1=0;
                                        if(res_CC3k[j]<=0.05) // 0 is the cut off value for structure in the out of density map
                                        {
                                            int kx0=j-1;
                                            if(kx0<0) kx0=0;
                                            initx0=kx0;
                                            for(int i=j;i<pnum;i++)
                                            {
                                                if(res_CC3k[i]>0.05) 
                                                {
                                                    float randcx=(2.0*randf0and1()-1.0);                                
                                                    if(randcx<0) break;
                                                }
                                                j++;
                                            }
                                            int kx1=j+1;
                                            if(kx1>pnum-1) kx1 = pnum-1;
                                            initx1=kx1;
                                        }
                                        int dis_01 = initx1-initx0+1;
                                        if(dis_01>=5) // 7 is the number which judgement how many residue should be as fragment
                                        {
                                            sssegment smg;
                                            smg.init = initx0;
                                            smg.term = initx1;
                                            segm.push_back(smg);
                                        }
                                    }
                                }                                                        
                    /*        for(int j=0;j<pnum;j++)
                            {
                                bb[j].indn = 3*j+0;
                                bb[j].indca = 3*j+1;
                                bb[j].indc = 3*j+2;
                                string resng = pointsrenm[j];
                                char resny = resng[0];
                                bb[j].resid = aminoid(resny);                                     
                            }     
                        //    calcsse2(bb, pnum, fin_maty);
                            calcssennhoc(bb,pnum,fin_maty);    */                      
                        //    rhoC0 =rhoC0;
                        //    cout<<" new_dE RRRRRRRRRRRR: "<<new_dE<<endl; 
                        }                
                        else
                        {
                        //    float tmpx=rand()/double(RAND_MAX);
                        //    cout<<"randx: "<<randtx<<endl;
                            float tmpx=randf0and1();
                            float mc_v = exp(-dE/(KT));
                            if(tmpx < mc_v)
                            {
                                int tmp_tmx=0;
                                for(int j=0;j<pnum;j++)
                                {
                                    pointsBx[3*j+0] = fin_maty[3*j+0];
                                    pointsBx[3*j+1] = fin_maty[3*j+1];
                                    pointsBx[3*j+2] = fin_maty[3*j+2];
                                //    bb[j]=bby[tmp_tmx];
                                    tmp_tmx = tmp_tmx +1;
                                }
                //                cout<<"KT,new_dE: "<<KT<<" "<<new_dE<<endl;
                                for(int j=0;j<pnum;j++)
                                {
                                    bb[j].indn = 3*j+0;
                                    bb[j].indca = 3*j+1;
                                    bb[j].indc = 3*j+2;  
                                    string resng = pointsrenm[j];
                                    char resny = resng[0];
                                    bb[j].resid = aminoid(resny);                                       
                                }     
                            //    calcsse2(bb, pnum, tmp_mat); 
                                calcssennhoc(bb,pnum,fin_maty);  

                                // calculate CC between the density map and ervery residue
                                vector<sssegment>().swap(segm);
                                if(calres)
                                {
                                    res_CC = vector<float>(pnum,2.0);
                                //    vector<float> res_CC3(pnum,0.0);
                                    nhg=0;
                                    for(int j=0;j<pnum;j++)
                                    {
                                        vector<poseCoord> sing_res(1);
                                        sing_res[0]= pointsBx[3*j+1];
                                        res_CC[j] = theDensityMap.matchposex(sing_res);
                                        res_CC3k[j]=res_CC[j];
                                    //    res_CC3[j]=res_CC3k[j];
                                //        cout<<"j,CC: "<<j<<" "<<res_CC[j]<<endl; 
                                    }
                                    
                                    
                                /*    for(int j=0;j<pnum;j++)
                                    {
                                        for(int jp=-1;jp<2;jp++)
                                        {
                                            int jjp = j+ jp;
                                            if(jjp<0) jjp=0;
                                            if(jjp>pnum-1) jjp=pnum-1;
                                            res_CC3[j] = res_CC3[j] + res_CC[jjp];
                                        }
                                        res_CC3k[j] = res_CC3[j]/3.0;
                                //        if(res_CC3k[j]<0.05) nhg=nhg+1;
                                    }     */
                                //    w1=w11*(float(pnum-nhg)/float(pnum))*(float(pnum-nhg)/float(pnum));
                                //    w6=w66*(float(pnum)/float(pnum-nhg))*(float(pnum)/float(pnum-nhg));            
                                    for(int j=0;j<pnum;j++)
                                    {
                                        int initx0=0;
                                        int initx1=0;
                                        if(res_CC3k[j]<=0.05) // 0 is the cut off value for structure in the out of density map
                                        {
                                            int kx0=j-1;
                                            if(kx0<0) kx0=0;
                                            initx0=kx0;
                                            for(int i=j;i<pnum;i++)
                                            {
                                                if(res_CC3k[i]>0.05) 
                                                {
                                                    float randcx=(2.0*randf0and1()-1.0);                                
                                                    if(randcx<0) break;
                                                }
                                                j++;
                                            }
                                            int kx1=j+1;
                                            if(kx1>pnum-1) kx1 = pnum-1;
                                            initx1=kx1;
                                        }
                                        int dis_01 = initx1-initx0+1;
                                        if(dis_01>=5) // 7 is the number which judgement how many residue should be as fragment
                                        {
                                            sssegment smg;
                                            smg.init = initx0;
                                            smg.term = initx1;
                                            segm.push_back(smg);
                                        }
                                    }
                                }                                 
                        /*        for(int j=0;j<pnum;j++)
                                {
                                    bb[j].indn = 3*j+0;
                                    bb[j].indca = 3*j+1;
                                    bb[j].indc = 3*j+2;   
                                    string resng = pointsrenm[j];
                                    char resny = resng[0];
                                    bb[j].resid = aminoid(resny);                                      
                                }     
                            //    calcsse2(bb, pnum, fin_maty);  
                                calcssennhoc(bb,pnum,fin_maty);        */                     
                            //    rhoC0 =rhoC0;
                            //    cout<<" XXX new_dE: "<<new_dE<<endl; 
                            }       
                        }                         
            }                     
        }
    }

/*    for(int j=0;j<pnum;j++)
    {
        pointsBx[3*j+0] = best_model[3*j+0];
        pointsBx[3*j+1] = best_model[3*j+1];
        pointsBx[3*j+2] = best_model[3*j+2];
    //    bb[j]=bby[tmp_tmx];
    }  */

    vector<poseCoord> pointsBy = pointsBx;
    cout<< "TTTTTT"<<endl;
    vector<poseCoord>().swap(best_model);
    best_model = pointsBx;
    best_E = 100000.0;
    int rec_num=50; // 5 
    int rec_numx = 30; // 10
    calres = true;
    w1=500.0;
    w6=60.0;
    w11=500.0;
    w66=60.0;
    w1=w11*(float(pnum-nhg)/float(pnum));
    w6=w66*(float(pnum)/float(pnum-nhg)); 
    for(int jjj=0;jjj<rec_num;jjj++)
    {
        float KT0s = 3.0;//;0.0001;
        float KT0e = 0.5;//;0.00001;
        float KT0x= pow(float(KT0e/KT0s),float(float(jjj)/float(rec_num)));
        float KT0xx = KT0s*float(KT0x);
        float KTns = 0.1;//;0.00005;
        float KTne = 0.01;//;0.0000001;
        float KTnx= pow(float(KTne/KTns),float(float(jjj)/float(rec_num)));
        float KTnxx = KTns*float(KTnx);   
        vector<float> Acprate(rec_numx,0.0);
        for(int jjb=0;jjb<rec_numx;jjb++)
        {
    //        theDensityMap.select_points( tmp_mat, theDensityMap.density);
    //        int nnum = theDensityMap.points_to_search_.size();            
            int Aprate0=0;
            int Aprate1=0;
            float KTx= pow(float(KTnxx/KT0xx),float(float(jjb)/float(rec_numx)));
            KT = KT0xx*float(KTx); 
        //    KT = 0.001;
        //    cout<<"KT: "<<KT<<endl;
        //    string elt_i;
            int allnum=500;
            int pnum_tmp=pnum;
            int dm_s=0; // start point for each domain
            int dm_e=pnum; // end point for each domain
            vector<sssegment>().swap(segm);
            res_CC = vector<float>(pnum,0.0);
            if(calres)
            {
                vector<float> res_CC3(pnum,0.0);
                nhg=0;
                for(int j=0;j<pnum;j++)
                {
                    vector<poseCoord> sing_res(1);
                    sing_res[0]= pointsBx[3*j+1];
                    res_CC[j] = theDensityMap.matchposex(sing_res);
                    if(res_CC[j]<0.05) nhg=nhg+1;
        //            cout<<"j,res_CC: "<<j<<" "<<res_CC[j]<<endl; 
                    vector<int>().swap(num_rand);
                    if(res_CC[j]<0.05)
                    {
                        num_rand.push_back(j);
                    }
                }    
                w1=w11*(float(pnum-nhg+1)/float(pnum));
                w6=w66*(float(pnum)/float(pnum-nhg+1));    
                for(int j=0;j<pnum;j++)
                {
                    int initx0=0;
                    int initx1=0;
                    if(res_CC[j]<0.05) // 0 is the cut off value for structure in the out of density map
                    {
                        int kx0=j;
                        if(kx0<0) kx0=0;
                        initx0=kx0;
                        for(int i=j;i<pnum;i++)
                        {
                            if(res_CC[i]>0.05) 
                            {
//                                float randcx=(2.0*randf0and1()-1.0);
//                                if(randcx<0.0) break;
                                break;
                            }
                            j++;
                        }
                        int kx1 = j;
                        if(kx1>pnum-1) kx1 = pnum-1;
                        initx1=kx1;
                    } 
                    int dis_01 = initx1-initx0+1;
                    if(dis_01>5) // 7 is the number which judgement how many residue should be as fragment
                    {
                        sssegment smg;
                        smg.init = initx0;
                        smg.term = initx1;
                        segm.push_back(smg);
                    }
                }                
            }          
                bool TF_segmx = false;
                if(segm.size()>0)
                {
                    TF_segmx=true;
                }                  
            for(int tt=0;tt<allnum;tt++)
            {       
        //        int randpx = 36;
        //        float KT0 =0.001;
        //        float KTn =0.0001;
                float sumis=0.0; // intersection
                float sumu=0.0; // union                
                vector<poseCoord > tmp_mat;
                vector<float> cartX2;
                vector<float> fracX2;        
                int rand_a2=0;
                int rand_b2=0;
        //        KT=0.001;
        //        float angle0 =180.0;
        //        float angle0 =2.0;
        //        float anglen =10;
        //        float anglen =1.0;
        //        float anglex= pow(float(anglen/angle0),float(float(tt)/float(allnum)));
        //        ang = angle0*float(anglex);
                ang=90.0;
        //        int randp = rand()%15+30;
                int recT = 0;   
                int rec_t = 0 ;  
                int rand_a = 0;
                int rand_b = 0; 
                int flagx=0 ;
            //    cout<<"hhhh"<<endl;
                vector<poseCoord>().swap(tmp_mat);
                //    int randtx = rand()%pnum_tmp+dm_s;
            //    int randtx = randIntCustom(0,pnum_tmp-1)+dm_s;
                int randtx =0;
                if(tt%5==0)
                {
                    int num_randss=num_rand.size();
                    if(num_randss>0)
                    {
                        int numrndst=randIntCustom(0,num_randss-1);
                        randtx=num_rand[numrndst];
                    } else
                    {
                        randtx = randIntCustom(0,pnum_tmp-1)+dm_s;
                    }
                    
                }else
                {
                    randtx = randIntCustom(0,pnum_tmp-1)+dm_s;
                }

                int randt=0;
                if(randtx>=(dm_e-6)||randtx<=dm_s+6)
                {
                    int jmx=randtx,jnx=randtx;
                    if(randtx>=(dm_e-6))
                    {
                        jmx=randtx;
                        jnx= dm_e-1;
                        rand_a=jmx;
                        rand_b=dm_e-1;                        
                    }else
                    {
                        jmx=0+dm_s;
                        jnx= randtx;  
                        rand_a=0+dm_s;
                        rand_b=jnx;                         
                    }
                    rand_a2 = rand_a;
                    rand_b2 = rand_b;                     
                }else
                {
                    randt = randIntCustom(0,4)+3;
                    int jmx=randtx,jnx=randtx;
                    float randt_num = 2.0*randf0and1()-1.0;
                    if(randt_num>0)
                    {
                        for(int j=0;j<randt;j++)
                        {
                            jnx=randtx+j;
                            if(jnx>=(dm_e-1)) break;
                        }
                        rand_a=jmx;
                        rand_b=jnx;                                                
                    }else
                    {
                        for(int j=0;j<randt;j++)
                        {
                            jmx=randtx-j;
                            if(jmx<=dm_s) break;
                        } 
                        rand_a=jmx;
                        rand_b=jnx;                                               
                    }
                    rand_a2 = rand_a-2;
                    rand_b2 = rand_b+2;                      
                }
            /*    if(randtx>=(dm_e-6)||randtx<=dm_s+6)
                {        
                    int jmx=randtx,jnx=randtx;            
                    if(randtx>=(dm_e-6))
                    {
                        jmx=randtx;
                        jnx= dm_e-1;
                        if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E')
                        {
                            for(int j=0;j<pnum_tmp;j++)
                            {
                                jmx=randtx-j;
                                if(jmx<=dm_s)
                                {
                                    jmx=dm_s;
                                    break;
                                }
                                if(numtoss(bb[jmx].sst)=='C') break;
                            }                 
                        }  
                        rand_a=jmx;
                        rand_b=dm_e-1;
                    } else
                    {
                        jmx=0+dm_s;
                        jnx= randtx;                
                        if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E')
                        {
                            for(int j=0;j<pnum_tmp;j++)
                            {
                                jnx=randtx+j;
                                if(jnx>=dm_e-1)
                                {
                                    jnx=dm_e-1;
                                    break;
                                }
                                if(numtoss(bb[jnx].sst)=='C') break;
                            }                  
                        }                     
                        rand_a=0+dm_s;
                        rand_b=jnx;                  
                    }
                    rand_a2 = rand_a;
                    rand_b2 = rand_b;               
                } else  
                {
                    if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E')
                    {
                    //    randt=rand()%4+6;
                        randt=randIntCustom(0,3)+6;
                        int jmx=randtx,jnx=randtx;
                        for(int j=0;j<randt;j++)
                        {
                            jnx=randtx+j;
                            if(jnx>=(dm_e-1)) break;
                            if(numtoss(bb[jnx].sst)=='C')
                            {
                                jnx=randtx+j-1;
                                break;
                            }
                        }
                        if((jnx-randtx+1)<randt)
                        {
                            for(int j=0;j<(randt-(jnx-randtx));j++)
                            {
                                jmx=randtx-j;
                                if(jmx<=dm_s) break;
                                if(numtoss(bb[jmx].sst)=='C') 
                                {
                                    jmx=randtx-j+1;
                                    break;
                                }
                            }
                        }
                        rand_a=jmx;
                        rand_b=jnx;
                    } else
                    {
                    //    randt = rand()%4+3;
                        randt = randIntCustom(0,4)+3;
                        int jmx=randtx,jnx=randtx;
                        for(int j=0;j<randt;j++)
                        {
                            jnx=randtx+j;
                            if(jnx>=(dm_e-1)) break;
                            if(numtoss(bb[jnx].sst)!='C') 
                            {
                                jnx=randtx+j-1;
                                break;
                            }
                        }
                        if((jnx-randtx+1)<randt)
                        {
                            for(int j=0;j<(randt-(jnx-randtx));j++)
                            {
                                jmx=randtx-j;
                                if(jmx<=dm_s) break;
                                if(numtoss(bb[jmx].sst)!='C') 
                                {
                                    jmx=randtx-j+1;
                                    break;
                                }
                            }
                        }
                        rand_a=jmx;
                        rand_b=jnx;                    
                    }
                    rand_a2 = rand_a-2;
                    rand_b2 = rand_b+2;              
                }         */
           /*     for(int i=0;i<segm.size();i++)
                {
                    int io1 = segm[i].init;
                    int io2 = segm[i].term;
                    if(rand_a<io2 && rand_a>io1) rand_a = io1;
                    if(rand_b<io2 && rand_b>io1) rand_b = io2;
                }            */

        /*        vector<int> segnum;
                vector<int>().swap(segnum);
                if(TF_segmx)
                {
             //       float randt_vs = 2.0*randf0and1()-1.0;
                //    if(randt_vs>=0)
                //    {
                        int flagtt=0;
                        for(int i=0;i<segm.size();i++)
                        {
                            flagtt=0;
                            int io1 = segm[i].init;
                            int io2 = segm[i].term;
                         //   if(rand_b<io2 && rand_a>io1) 
                         //   {
                         //       rand_a = io1;
                         //       rand_b = io2;
                         //   }                           
                            if(rand_a<io2 && rand_a>io1) {rand_a = io1;flagtt=1;}
                            if(rand_b<io2 && rand_b>io1) {rand_b = io2;flagtt=1;}
                            if(flagtt==1){
                                segnum.push_back(i);
                            }                            
                        }
                //    } 
                }         */                
            //    rand_a2 = rand_a-2;
            //    rand_b2 = rand_b+2;
                if(rand_a>rand_b)
                {
                    int tmp_rand = rand_b;
                    rand_b = rand_a;
                    rand_a = tmp_rand;
                }  
                if(rand_a <0) rand_a = 0;
                if(rand_b >(pnum-1)) rand_b= pnum-1;

                old_dE = 0.0;
        //        if(rand_b2>(pnum-1)||rand_a2<0) continue;
                if(rand_b2>(pnum-1)) rand_b2 = (pnum-1);
                if(rand_a2<0) rand_a2 =0;
                //    cout<<"dddd"<<endl; 
                if(abs(rand_b-rand_a+1)<3) continue;
                 //        cout<<"eee"<<endl;
    
        //        cout<<"rand a b: "<<rand_a<<" "<<rand_b<<endl;          
        //                if( abs(rand_a -rand_b+1)<4) continue;
                for(int j=0;j<pnum;j++)
                {
                    tmp_mat.push_back(pointsBx[3*j+0]);
                    tmp_mat.push_back(pointsBx[3*j+1]);
                    tmp_mat.push_back(pointsBx[3*j+2]);
                }               
                float old_Ehbond=0.0;

        //        if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E')
                int bn=1;
                if(bn==1)
                {
            /*        for(int j=0;j<pnum;j++)
                    {
                        bb[j].indn = 3*j+0;
                        bb[j].indca = 3*j+1;
                        bb[j].indc = 3*j+2;
                        string resng = pointsrenm[j];
                        char resny = resng[0];
                        bb[j].resid = aminoid(resny);                             
                    }     
                //    calcsse2(bb, pnum, pointsBx);  
                    calcssennhoc(bb,pnum,pointsBx);           */   
                    int numseq=pnum;
                    vector<point3f> decstr(numseq);                
                    tmp_tt=0;
                    for(int j=0;j<pnum;j++)
                    {
                        decstr[tmp_tt].ss2 = numtoss(bb[j].sst);// predicted secondary structure
                        decstr[tmp_tt].ssm = numtoss(bb[j].sst);
                        decstr[tmp_tt].x = pointsBx[3*tmp_tt+1].x_[0];
                        decstr[tmp_tt].y = pointsBx[3*tmp_tt+1].x_[1];
                        decstr[tmp_tt].z = pointsBx[3*tmp_tt+1].x_[2];
                        decstr[tmp_tt].ptn.x = pointsBx[3*tmp_tt+0].x_[0];
                        decstr[tmp_tt].ptn.y = pointsBx[3*tmp_tt+0].x_[1];
                        decstr[tmp_tt].ptn.z = pointsBx[3*tmp_tt+0].x_[2];
                        decstr[tmp_tt].ptc.x = pointsBx[3*tmp_tt+2].x_[0];
                        decstr[tmp_tt].ptc.y = pointsBx[3*tmp_tt+2].x_[1];
                        decstr[tmp_tt].ptc.z = pointsBx[3*tmp_tt+2].x_[2];
                        tmp_tt=tmp_tt+1;
                    } 
                //    old_Ehbond=energyhbondcanc(decstr,numseq);
                    old_Ehbond=energyhbondnhoc2(decstr,numseq);
                    vector<point3f>().swap(decstr);  
                }            
        /*            float old_Ehbond=0.0;
                    int numseq=abs(rand_b-rand_a+1);
                    if(numseq<3) continue;
                    vector<point3f> decstr(numseq);

                    if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E')
                    {
                        tmp_tt=0;
                        for(int j=rand_a;j<=rand_b;j++)
                        {
                            decstr[tmp_tt].ss2=numtoss(bb[j].sst);
                            decstr[tmp_tt].x=tmp_mat[3*tmp_tt+1].x_[0];
                            decstr[tmp_tt].y=tmp_mat[3*tmp_tt+1].x_[1];
                            decstr[tmp_tt].z=tmp_mat[3*tmp_tt+1].x_[2];
                            decstr[tmp_tt].ptn.x=tmp_mat[3*tmp_tt+0].x_[0];
                            decstr[tmp_tt].ptn.y=tmp_mat[3*tmp_tt+0].x_[1];
                            decstr[tmp_tt].ptn.z=tmp_mat[3*tmp_tt+0].x_[2];
                            decstr[tmp_tt].ptc.x=tmp_mat[3*tmp_tt+2].x_[0];
                            decstr[tmp_tt].ptc.y=tmp_mat[3*tmp_tt+2].x_[1];
                            decstr[tmp_tt].ptc.z=tmp_mat[3*tmp_tt+2].x_[2];
                            tmp_tt=tmp_tt+1;
                        } 
                        old_Ehbond=energyhbondcanc(decstr,numseq);
                    }  */
                //    cout<<"old_Ehbond: "<<old_Ehbond<<" ";

                    // ss information
            /*        float old_ss_E=0.0;
                    int CA_num=0;
                    vector<vector<float> > old_CA_cor;
                    int tm_i1=0;               
                    for(int j=0;j<tmp_mat.size();j++)
                    { 
                        if(tmp_mat[j].elt_ == "CA")
                        {
                            old_CA_cor.push_back(tmp_mat[j].x_);
                            CA_num = CA_num+1;
                            tm_i1=tm_i1+1; 
                        }                  
                    }
                    vector<string > old_ss_inf(CA_num);
                    make_secx(old_CA_cor,CA_num,old_ss_inf);
                    int ss_ab_num=0; // helix and b strand number
                    int ss_c_num =0; // coil number
        //                float old_ss_ab=0.0;
        //                float old_ss_c=0.0;                
                    for(int j=0;j<old_ss_inf.size();j++)
                    {
                        if(old_ss_inf[j]=='H')
                        {
                            old_ss_E = old_ss_E +0.1;
        //                        old_ss_ab = old_ss_ab + 0.4;
        //                        ss_ab_num = ss_ab_num +1 ;
                        }
                        else if(old_ss_inf[j]=='E')
                        {
                            old_ss_E = old_ss_E +0.1;
        //                        old_ss_ab = old_ss_ab + 0.4;
        //                        ss_ab_num = ss_ab_num +1;
                        }
                        else{
                            old_ss_E = old_ss_E +0.2;
        //                        old_ss_c = old_ss_c + 0.5;
        //                        ss_c_num = ss_c_num +1;
                        }                    
                    }
        //                old_ss_E = (old_ss_ab*ss_ab_num+old_ss_c*ss_c_num)/(ss_ab_num+ss_c_num);
                    old_ss_E = old_ss_E/(CA_num); */

                    ObjexxFCL::FArray3D< float > rhoC2;
                    ObjexxFCL::FArray3D< float > inv_rho_mask2;
                    rhoC2.dimension(Denx , Deny , Denz);
                    inv_rho_mask2.dimension(Denx , Deny , Denz);
                    for ( int t=0; t<Denxyz; ++t ) {
                        rhoC2[t]=0.0;
                        inv_rho_mask2[t]=1.0;
                    } 
                    vector<float> del_ij(3,0.0);
                    vector<float> atm_j(3,0.0); 
                    int tm_i1=0;          
                    del_ij=vector<float>(3,0.0);
                    atm_j=vector<float>(3,0.0);
                    tm_i1 =0;
                    vector<float> atm_idy(3,0.0);
                    bond0max=0.0;bond1max=0.0;bond2max=0.0;
                    bond0min=10000.0;bond1min=10000.0;bond2min=10000.0;                
            //            atm_idy=vector<vector<float> > ((rand_b-rand_a+1),vector<float>(3,0.0));
                    for(int j=3*rand_a2;j<=3*rand_b2+2;j++)
                    {
                        if(tmp_mat[j].elt_ == "CA")
                        {
                            vector<float> cartX1;
                            vector<float> fracX1; 
                        //    elt_i = tmp_mat[3*tm_i1+1].elt_;
                            elt_i = tmp_mat[j].elt_;
                            elt_i = elt_i[0];
                            OneGaussianScattering sig_j = get_A( elt_i );
                            k = sig_j.k( theDensityMap.effectiveB );
                            C = sig_j.C( k );
                            if ( C < 1e-6 ) continue; 

                    //        cartX1 = tmp_mat[3*tm_i1+1].x_;
                            cartX1 = tmp_mat[j].x_;
                    //        tm_i1 = tm_i1 + 1;
                            MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);
                            atm_idy[0] = pos_mod (double(fracX1[0]*grid[0] - origin[0] + 1) , (double)grid[0]);
                            atm_idy[1] = pos_mod (double(fracX1[1]*grid[1] - origin[1] + 1) , (double)grid[1]);
                            atm_idy[2] = pos_mod (double(fracX1[2]*grid[2] - origin[2] + 1) , (double)grid[2]);                 
                            for(int z=1;z<=Denz;z++)
                            {
                                atm_j[2] = z;
                                del_ij[2] =(atm_idy[2]-atm_j[2])/grid[2];
                                if ( del_ij[2] > 0.5 ) del_ij[2]-=1.0;
                                if ( del_ij[2] < -0.5 ) del_ij[2]+=1.0;                    
                                del_ij[0] = del_ij[1] = 0.0;
                                vector<float> frac_tmpz;
                                MatrixTimesTransVector(theDensityMap.f2c,del_ij,frac_tmpz);
                            //    cout<<"ATOM: "<<theDensityMap.ATOM_MASK_PADDING+theDensityMap.ATOM_MASK<<endl;
                                if(square_len(frac_tmpz)> (padding+atom_m)*(padding+atom_m) ) continue; 
                                if(z<bond2min) bond2min = z;
                                if(z>bond2max) bond2max = z;              
                                for(int y=1;y<=Deny;y++)
                                {
                                    atm_j[1] = y;
                                    del_ij[1] = (atm_idy[1] - atm_j[1])/grid[1];
                                    // wrap-around??
                                    if ( del_ij[1] > 0.5 ) del_ij[1]-=1.0;
                                    if ( del_ij[1] < -0.5 ) del_ij[1]+=1.0;                        
                                    del_ij[0] = 0.0;
                                    vector<float> frac_tmpy;
                                    MatrixTimesTransVector(theDensityMap.f2c,del_ij,frac_tmpy);
                                    if(square_len(frac_tmpy)> (padding+atom_m)*(padding+atom_m) ) continue;
                                    if(y>bond1max) bond1max = y;
                                    if(y<bond1min) bond1min = y;                      
                                    for(int x=1;x<=Denx;x++)
                                    {
                                        atm_j[0] = x;
                                        del_ij[0] = (atm_idy[0] - atm_j[0])/grid[0];
                                        // wrap-around??
                                        if ( del_ij[0] > 0.5 ) del_ij[0]-=1.0;
                                        if ( del_ij[0] < -0.5 ) del_ij[0]+=1.0;                            
                                        vector<float> cart_del_ij2;
                                        MatrixTimesTransVector(theDensityMap.f2c,del_ij,cart_del_ij2);
                                        float d2 = square_len(cart_del_ij2);
                                        if(d2 > (padding+atom_m)*(padding+atom_m) ) continue;
                                    
                                        float atm = C*exp(-k*d2);
                                        float sigmoid_msk = exp( d2 - (atom_m)*(atom_m)  );
                                    //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK));
                                    //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK));
                                        float inv_msk = 1.0/(1.0+sigmoid_msk);
                        //                rhoC2(x,y,z) += atm;
                        //                rhoC2(x,y,z) += atm;
                                        rhoC2(x,y,z) += atm;
                         //               if(rhoC01(x,y,z)<0) rhoC01(x,y,z)=0;
                                        inv_rho_mask2(x,y,z) *= (1.0 - inv_msk);

                                        if(x>bond0max) bond0max = x;  
                                        if(x<bond0min) bond0min = x;                                                     

            //                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
            //                                mask2(x,y,z) = 1.0; // problem?
            //                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
            //                                    float atm = C*exp(-k*d2);
            //                                    rhoC(x,y,z) += atm;
            //                                }
            //                            }
                                    }
                                }
                            }
                            tm_i1 = tm_i1 + 1;
                        }
                    }

         //           for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
         //               rhoC2[x] = rhoC0[x];
        //            }            
        //            float max_rhc = 0.0;
        //            for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
        //                if(rhoC2[x]>max_rhc) max_rhc = rhoC2[x] ;
        //            }
            /*        for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
                //        cout<<theDensityMap.density[x]<<" ";
                        if(rhoC2[x]<(3.0/5.0)*max_rhc) rhoC2[x]=0; 
                    }    */
        //            rhoC2 = rhoC0;         
                    float sumC_i2=0.0, sumO_i2=0.0, sumCO_i2=0.0, vol_i2=0.0, CC_i2=0.0;
                    float sumO2_i2=0.0, sumC2_i2=0.0, varC_i2=0.0, varO_i2=0.0;
                    float clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
            /*        for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
        //                tmp_den = 0.0; 
        //                tmp_rhc=0.0;
          //              if(theDensityMap.density[x]>(2.0/3.0)*max_den) tmp_den=theDensityMap.density[x];
        //                if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
        //                if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                        // fetch this point
                        clc_x2 = rhoC2[x];
        //                clc_x2 = tmp_rhc;
                        obs_x2 = theDensityMap.density[x];
        //                obs_x2 = tmp_den;
                        eps_x2 = 1-inv_rho_mask2[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                    //    eps_x2 = 1.0;

                        // SMOOTHED
                        sumCO_i2 += eps_x2*clc_x2*obs_x2;
                        sumO_i2  += eps_x2*obs_x2;
                        sumO2_i2 += eps_x2*obs_x2*obs_x2;
                        sumC_i2  += eps_x2*clc_x2;
                        sumC2_i2 += eps_x2*clc_x2*clc_x2;
                        vol_i2   += eps_x2;
                    }
                    varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
                    varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
                    if ( varC_i2 == 0 || varO_i2 == 0 ) {
                        CC_i2 = 0;
                    } else {
                        CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
                    }  */
                    for(int x=int(bond0min);x<=bond0max;x++)
                    {
                       for(int y=int(bond1min);y<=bond1max;y++)
                        {
                            for(int z=int(bond2min);z<=bond2max;z++)
                            {
                //                tmp_den = 0.0; 
                //                tmp_rhc=0.0;
                  //              if(theDensityMap.density[x]>(2.0/3.0)*max_den) tmp_den=theDensityMap.density[x];
                //                if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
                //                if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                                // fetch this point
                                clc_x2 = rhoC2(x,y,z);
                //                clc_x2 = tmp_rhc;
                                obs_x2 = theDensityMap.density(x,y,z);
                //                obs_x2 = tmp_den;
                                eps_x2 = 1-inv_rho_mask2(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                            //    eps_x2 = 1.0;

                                // SMOOTHED
                                sumCO_i2 += eps_x2*clc_x2*obs_x2;
                                sumO_i2  += eps_x2*obs_x2;
                                sumO2_i2 += eps_x2*obs_x2*obs_x2;
                                sumC_i2  += eps_x2*clc_x2;
                                sumC2_i2 += eps_x2*clc_x2*clc_x2;
                                vol_i2   += eps_x2;
                            }
                        }
                    }          
                    varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
                    varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
                    if ( varC_i2 == 0 || varO_i2 == 0 ) {
                        CC_i2 = 0;
                    } else {
                        CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
                    }           
            /*        sumis=0.0;
                    sumu=0.0;
                    for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
                        tmp_den = 0.0; 
                        tmp_rhc=0.0;
                        if(theDensityMap.density[x]>(3.0/5.0)*max_den) tmp_den=1.0;
                        if(rhoC2[x]>(3.0/5.0)*max_rhc) tmp_rhc=1.0;
        //                if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                        // fetch this point
                        sumis=sumis+tmp_den*tmp_rhc;
                        sumu=sumu+tmp_den+tmp_rhc-tmp_den*tmp_rhc;
                    }
                    if ( sumu == 0 ) {
                        CC_i2 = 0;
                    } else {
                        CC_i2 = sumis/sumu;
                    }            */
        //            vector<vector<float> > bondx(3,vector<float>(2,0.0));
           /*         float bond0max=0.0,bond1max=0.0,bond2max=0.0;
                    float bond0min=1000.0,bond1min=1000.0,bond2min=1000.0;
        //            vector<vector<float>> atm_idy((rand_b-rand_a+1),vector<float>(3,0.0));
                    atm_idy=vector<vector<float>>((rand_b-rand_a+1),vector<float>(3,0.0));
                    tm_i1=0; 
                    for(int j=0;j<tmp_mat.size();j++)
                    {
         //               tmp_mat.push_back(pointsBx[j]);
                        if(tmp_mat[j].elt_ == "CA")
                        {
                            vector<float> cartX1;
                            vector<float> fracX1;
            //                cartX1 = tmp_mat[3*tm_i1+1].x_;
                            cartX1 = tmp_mat[j].x_;
                            MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);
                            atm_idy[tm_i1][0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                            atm_idy[tm_i1][1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                            atm_idy[tm_i1][2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);
                            if(atm_idy[tm_i1][0]>bond0max) bond0max = atm_idy[tm_i1][0];
                            if(atm_idy[tm_i1][1]>bond1max) bond1max = atm_idy[tm_i1][1];
                            if(atm_idy[tm_i1][2]>bond2max) bond2max = atm_idy[tm_i1][2];
                            if(atm_idy[tm_i1][0]<bond0min) bond0min = atm_idy[tm_i1][0];
                            if(atm_idy[tm_i1][1]<bond1min) bond1min = atm_idy[tm_i1][1];
                            if(atm_idy[tm_i1][2]<bond2min) bond2min = atm_idy[tm_i1][2];
                            tm_i1 = tm_i1 + 1;
                        }
                    }  */
            /*        sumis=0.0;
                    sumu=0.0;            
                    for(int x=int(bond0min);x<=bond0max;x++)
                    {
                       for(int y=int(bond1min);y<=bond1max;y++)
                        {
                            for(int z=int(bond2min);z<=bond2max;z++)
                            {
                                tmp_den = 0.0; 
                                tmp_rhc=0.0;
                                if(theDensityMap.density(x,y,z)>=(3.0/5.0)*max_den) tmp_den=1.0;
                                if(rhoC2(x,y,z)>(3.0/5.0)*max_rhc) tmp_rhc=1.0;
                //                if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                                // fetch this point
                                sumis=sumis+tmp_den*tmp_rhc;
                                sumu=sumu+tmp_den+tmp_rhc-tmp_den*tmp_rhc;
                            }
                        }
                    }          
                    if ( sumu == 0 ) {
                        CC_i2 = 0;
                    } else {
                        CC_i2 = sumis/sumu;
                    }  */

        /*            float bond0max=0.0,bond1max=0.0,bond2max=0.0;
                    float bond0min=1000.0,bond1min=1000.0,bond2min=1000.0;
                    vector<vector<float>> atm_idy((rand_b-rand_a+1),vector<float>(3,0.0));
                    tm_i1=0; 
                    for(int j=0;j<tmp_mat.size();j++)
                    {
         //               tmp_mat.push_back(pointsBx[j]);
                        if(tmp_mat[j].elt_ == "CA")
                        {
                            vector<float> cartX1;
                            vector<float> fracX1;
            //                cartX1 = tmp_mat[3*tm_i1+1].x_;
                            cartX1 = tmp_mat[j].x_;
                            MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);
                            atm_idy[tm_i1][0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                            atm_idy[tm_i1][1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                            atm_idy[tm_i1][2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);
                            if(atm_idy[tm_i1][0]>bond0max) bond0max = atm_idy[tm_i1][0];
                            if(atm_idy[tm_i1][1]>bond1max) bond1max = atm_idy[tm_i1][1];
                            if(atm_idy[tm_i1][2]>bond2max) bond2max = atm_idy[tm_i1][2];
                            if(atm_idy[tm_i1][0]<bond0min) bond0min = atm_idy[tm_i1][0];
                            if(atm_idy[tm_i1][1]<bond1min) bond1min = atm_idy[tm_i1][1];
                            if(atm_idy[tm_i1][2]<bond2min) bond2min = atm_idy[tm_i1][2];
                            tm_i1 = tm_i1 + 1;
                        }
                    } 
                    sumis=0.0;
                    sumu=0.0;            
                    for(int x=int(bond0min);x<=bond0max;x++)
                    {
                       for(int y=int(bond1min);y<=bond1max;y++)
                        {
                            for(int z=int(bond2min);z<=bond2max;z++)
                            {
                                clc_x2 = rhoC2(x,y,z);
                //                clc_x2 = tmp_rhc;
                                obs_x2 = theDensityMap.density(x,y,z);
                //                obs_x2 = tmp_den;
                //                eps_x2 = 1-inv_rho_mask2[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                                eps_x2 = 1.0;

                                // SMOOTHED
                                sumCO_i2 += eps_x2*clc_x2*obs_x2;
                                sumO_i2  += eps_x2*obs_x2;
                                sumO2_i2 += eps_x2*obs_x2*obs_x2;
                                sumC_i2  += eps_x2*clc_x2;
                                sumC2_i2 += eps_x2*clc_x2*clc_x2;
                                vol_i2   += eps_x2;
                            }
                        }
                    }          
                    varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
                    varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
                    if ( varC_i2 == 0 || varO_i2 == 0 ) {
                        CC_i2 = 0;
                    } else {
                        CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
                    }  */

        /*                for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
                            clc_x2 = rhoC2[x];
                            obs_x2 = theDensityMap.density[x];
                            eps_x2 = 1-inv_rho_mask2[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal

                            // SMOOTHED
                            sumCO_i2 += eps_x2*clc_x2*obs_x2;
                            sumO_i2  += eps_x2*obs_x2;
                            sumO2_i2 += eps_x2*obs_x2*obs_x2;
                            sumC_i2  += eps_x2*clc_x2;
                            sumC2_i2 += eps_x2*clc_x2*clc_x2;
                            vol_i2   += eps_x2;
                        }
                        sumO_i2 = sumO_i2/vol_i2;
                        sumC_i2 = sumC_i2/vol_i2;
                        sumCO_i2 = sumCO_i2/vol_i2;
                        sumC2_i2 = sumC2_i2/vol_i2;
                        sumO2_i2 = sumO2_i2/vol_i2;
                        varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2);
                        varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2) ;
                        if ( varC_i2 == 0 || varO_i2 == 0 ) {
                            CC_i2 = 0;
                        } else {
                            CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2) / sqrt( varC_i2 * varO_i2 );
                        } */
                /*        for(int x=bond0min;x<=bond0max;x++)
                        {
                           for(int y=bond1min;y<=bond1max;y++)
                            {
                                for(int z=bond2min;z<=bond2max;z++)
                                {
                                    clc_x2 = rhoC2(x,y,z);
                    //                clc_x2 = tmp_rhc;
                                    obs_x2 = theDensityMap.density(x,y,z);
                    //                obs_x2 = tmp_den;
                    //                eps_x2 = 1-inv_rho_mask2[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                                    eps_x2 = 1.0;

                                    // SMOOTHED
                                    vol_i2 += eps_x2;
                                    sumC_i2 += (eps_x2*clc_x2*clc_x2);
                                    sumO_i2 += (eps_x2*obs_x2*obs_x2);
                                    sumCO_i2 += (eps_x2*clc_x2*obs_x2);
                                }
                            }
                        }          
                        if ( sumC_i2 == 0 || sumO_i2 == 0 ) {
                            CC_i2 = 0;
                        } else {
                            CC_i2 = sumCO_i2/(sqrt(sumC_i2+sumO_i2)); ;
                        }             */
            /*            for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x )
                        {
                            eps_x2 = 1-inv_rho_mask2[x];
                            vol_i2   += eps_x2;
                            sumC_i2 += (eps_x2*rhoC2[x]*rhoC2[x]);
                            sumO_i2 += (eps_x2*theDensityMap.density[x]*theDensityMap.density[x]);
                            sumCO_i2 += (eps_x2*rhoC2[x]*theDensityMap.density[x]);
                        }
                //        sumCO_i2 = sumCO_i2/vol_i2;
                //        sumO_i2 = sumO_i2/vol_i2;
                //        sumC_i2 = sumC_i2/vol_i2;
                        CC_i2 = sumCO_i2/(sqrt(sumC_i2*sumO_i2));  */
                        old_CC =1.000- CC_i2; 
                //        cout<< "old_CC: "<<old_CC<<" ";

                        float dis_p=0.0;
                        bool tp_clashx= false;
                        old_vwd =0.0;
                        old_clash = 0.0;
                        for(int jj=3*rand_a2;jj<=3*rand_b2+2;jj++)
                        {
                            if(tmp_mat[jj].elt_=="CA")
                            {
        //                        float VR1= GetVdwRadius(tmp_mat[jj].elt_);
                                string VR1 = tmp_mat[jj].elt_ ;
                            //    for(int t=0;t<pointsBx.size();t++)
                                for(int t=3*dm_s+0;t<3*dm_e;t++)
                                {
                                    if(tmp_mat[t].elt_=="CA")
                                    {
                        //                if((jj+rand_a)==t) continue;
                                //        if(abs(jj+3*rand_a-t)<1) continue;
                                        if(abs(jj-t)<3) continue;
                               //         if(t>3*rand_a && t<(3*rand_a+tmp_mat.size())) continue;
        //                                float VR2= GetVdwRadius(pointsBx[t].elt_);
                                        string VR2 = tmp_mat[t].elt_ ;
                                        dis_p = Distance_point(tmp_mat[jj].x_,tmp_mat[t].x_);
        //                                if(dis_p < (VR1+VR2))
                            //            if(dis_p < 3.70)
                            //            {
                         //                    old_clash = old_clash + 1.0/ sqrt(dis_p+0.0001) ;
                        //                    old_clash = old_clash + exp( 3.75-dis_p) ;
                                            old_clash = old_clash + GetVdwEgC(VR1,VR2,dis_p);
                              //          }
                        //                old_vwd=old_vwd+GetVdwEg(VR1,VR2,dis_p);
                        //                old_vwd=old_vwd+GetVdwEgC(VR1,VR2,dis_p);
                    //                    cout<<"jj,t: "<<jj+rand_a<<" "<<t<<" dis: "<<VR1+VR2<<" "<<dis_p<<endl;
                    //                    cout<<"vwdeg: "<<GetVdwEg(VR1,VR2,dis_p)<<endl;;
                    //                    if(((dis_p)<2.0))
                    //                    {
                    //                        tp_clashx = true;
                    //                        break;
                    //                    } 
                    //                    if(tp_clashx)
                    //                    {
                    //                        cout<<"jj,t: "<<tmp_mat[jj].elt_<<" "<<pointsBx[t].elt_<<" dis: "<<VR1+VR2<<" "<<dis_p<<endl;
                    //                        tp_clashx=false;
                    //                    }
                    //                    cout<< "VR1+VR2,dis_p: "<<VR1+VR2<<" "<<dis_p<<endl;
                                    }
                                }
                            }
            //                if(tp_clashx) break;
                        }  

            /*        vector<point3f> decstrbondang(pnum);
                    for(int j=0;j<pnum;j++)
                    {
                        decstrbondang[j].ss2=numtoss(bb[j].sst);
                        decstrbondang[j].x=pointsBx[3*j+1].x_[0];
                        decstrbondang[j].y=pointsBx[3*j+1].x_[1];
                        decstrbondang[j].z=pointsBx[3*j+1].x_[2];
                        decstrbondang[j].ptn.x=pointsBx[3*j+0].x_[0];
                        decstrbondang[j].ptn.y=pointsBx[3*j+0].x_[1];
                        decstrbondang[j].ptn.z=pointsBx[3*j+0].x_[2];
                        decstrbondang[j].ptc.x=pointsBx[3*j+2].x_[0];
                        decstrbondang[j].ptc.y=pointsBx[3*j+2].x_[1];
                        decstrbondang[j].ptc.z=pointsBx[3*j+2].x_[2];
                        string resng = pointsrenm[j];
                        char resny = resng[0];
                        decstrbondang[j].iaa = aminoid(resny);
        //                decstr[i].iaa=pp.aminoid(decstr[i].aaa);                
        //                tmp_tt=tmp_tt+1;
                    }     
                    for(int j=1;j<pnum;j++)
                    {        
                        if(j==0)
                        {
                            decstrbondang[j].ang[0]=0.0;
                            decstrbondang[j].ang[1]=0.0;
                            decstrbondang[j].ang[2]=0.0;
                        }                   
                        point3d p12 = setv(decstrbondang[j-1].x-decstrbondang[j-1].ptc.x,decstrbondang[j-1].y-decstrbondang[j-1].ptc.y,decstrbondang[j-1].z-decstrbondang[j-1].ptc.z);
                        point3d p23 = setv(decstrbondang[j].ptn.x-decstrbondang[j-1].ptc.x,decstrbondang[j].ptn.y-decstrbondang[j-1].ptc.y,decstrbondang[j].ptn.z-decstrbondang[j-1].ptc.z);
                        decstrbondang[j].ang[0]= angv(p12,p23)*degrad;
                        p12 = setv(decstrbondang[j-1].ptc.x-decstrbondang[j].ptn.x,decstrbondang[j-1].ptc.y-decstrbondang[j].ptn.y,decstrbondang[j-1].ptc.z-decstrbondang[j].ptn.z);
                        p23 = setv(decstrbondang[j].x-decstrbondang[j].ptn.x,decstrbondang[j].y-decstrbondang[j].ptn.y,decstrbondang[j].z-decstrbondang[j].ptn.z);
                        decstrbondang[j].ang[1]= angv(p12,p23)*degrad;
                        p12 = setv(decstrbondang[j].ptn.x-decstrbondang[j].x,decstrbondang[j].ptn.y-decstrbondang[j].y,decstrbondang[j].ptn.z-decstrbondang[j].z);
                        p23 = setv(decstrbondang[j].ptc.x-decstrbondang[j].x,decstrbondang[j].ptc.y-decstrbondang[j].y,decstrbondang[j].ptc.z-decstrbondang[j].z);
                        decstrbondang[j].ang[2]= angv(p12,p23)*degrad;
        //                tmp_tt=tmp_tt+1;
                    }
            //        cout<<"FFF"<<endl; 
                    float old_bondangenergy =0.0;
                    Eenergybondangle(decstrbondang,pnum,old_bondangenergy);
                    cout<< "old_bondangenergy: "<< old_bondangenergy/pnum<<endl;                 */

        /*                vector<point3f> decstrbondang(numseq);
                        for(int j=0;j<pnum;j++)
                        {
                            decstrbondang[j].iaa = ;
                            decstrbondang[j].ang[0]= ;
                            decstrbondang[j].ss2=numtoss(bb[j].sst);
                            decstrbondang[j].x=pointsBx[3*j+1].x_[0];
                            decstrbondang[j].y=pointsBx[3*j+1].x_[1];
                            decstrbondang[j].z=pointsBx[3*j+1].x_[2];
                            decstrbondang[j].ptn.x=pointsBx[3*j+0].x_[0];
                            decstrbondang[j].ptn.y=pointsBx[3*j+0].x_[1];
                            decstrbondang[j].ptn.z=pointsBx[3*j+0].x_[2];
                            decstrbondang[j].ptc.x=pointsBx[3*j+2].x_[0];
                            decstrbondang[j].ptc.y=pointsBx[3*j+2].x_[1];
                            decstrbondang[j].ptc.z=pointsBx[3*j+2].x_[2];
            //                tmp_tt=tmp_tt+1;
                        }  */
            //            cout<<"tp_clashx: " <<tp_clashx<<endl;
            //            if(tp_clashx) continue;  
            //            old_dE = old_CC + old_vwd;
            //            old_dE = 30*old_CC + old_vwd/((pnum-1)*(rand_b-rand_a+1));
                   //     old_dE = 0.8*old_CC + 0.5*old_clash/abs(rand_b-rand_a+1)+0.5*old_Ehbond/abs(rand_b-rand_a+1);
            //            old_dE = 50.0*old_CC + 50.0*old_clash/abs(rand_b-rand_a+1)+50.0*old_Ehbond/abs(rand_b-rand_a+1);
                   //     old_dE = old_CC + old_clash/abs(rand_b2-rand_a2+1)+old_Ehbond/abs(pnum)/(1.0-(old_Ehbond/abs(pnum)));
                   //     old_dE = 500.0*old_CC + 100.0*old_clash/abs(rand_b2-rand_a2+1) + 100.0*(old_Ehbond/abs(pnum))/(1.0-(old_Ehbond/abs(pnum)));// + old_bondangenergy/pnum;
                        float Eres0=0.0;
                    /*    for(int j=0;j<pnum;j++)
                        {
                            Eres0 = Eres0 + Distance_point(pointsBx[3*j+1].x_,pointsBy[3*j+1].x_);
                        } */

                        for(int j=3*rand_a2;j<=3*rand_b2+2;j++)
                        {
                            for(int js=j-3;js<=j+3;js++)
                            {
                                if(js==j) continue;
                                if(js<0) continue;
                                if(js>(pnum-1)) continue;
                                float dx0 = Distance_point(tmp_mat[j].x_,tmp_mat[js].x_);
                                float dx1 = Distance_point(pointsBx0[j].x_,pointsBx0[js].x_);
                                Eres0 = Eres0 + (dx0-dx1)*(dx0-dx1);
                            }
                        }                        
                        Eres0 = Eres0/float(rand_b2-rand_a2+1); 

                        vector<point3f> decstrbondlen(pnum);
                        for(int j=0;j<pnum;j++)
                        {
                            decstrbondlen[j].ss2=numtoss(bb[j].sst);
                            decstrbondlen[j].ssm=numtoss(bb[j].sst);
                            decstrbondlen[j].stype=numtoss(bb[j].sst);
                            decstrbondlen[j].x=pointsBx[3*j+1].x_[0];
                            decstrbondlen[j].y=pointsBx[3*j+1].x_[1];
                            decstrbondlen[j].z=pointsBx[3*j+1].x_[2];
                            decstrbondlen[j].ptn.x=pointsBx[3*j+0].x_[0];
                            decstrbondlen[j].ptn.y=pointsBx[3*j+0].x_[1];
                            decstrbondlen[j].ptn.z=pointsBx[3*j+0].x_[2];
                            decstrbondlen[j].ptc.x=pointsBx[3*j+2].x_[0];
                            decstrbondlen[j].ptc.y=pointsBx[3*j+2].x_[1];
                            decstrbondlen[j].ptc.z=pointsBx[3*j+2].x_[2];
                            string resng = pointsrenm[j];
                            char resny = resng[0];
                            decstrbondlen[j].iaa = aminoid(resny);
            //                decstr[i].iaa=pp.aminoid(decstr[i].aaa);                
            //                tmp_tt=tmp_tt+1;
                        } 
                        str2tor(decstrbondlen,pnum,3);
                    //    str2torp(decstrbondlen,pnum,0,pnum-1);  
                /*        float old_fang=0.0;
                        float old_Ebondang=0.0;
                    //    cout<<"jjk: "<<endl;
                        old_Ebondang =(float) Eenergybondangle(decstrbondlen,pnum,old_fang);
                    //    cout<<"jjk1: "<<endl;
                        old_Ebondang = 0.30*old_Ebondang + old_fang;         */                                
                        float old_fene=0.0;
                        float old_Ebondlen=0.0;
                        old_Ebondlen = (float) energybondlength(decstrbondlen,pnum,old_fene);
                        old_Ebondlen = 0.5*old_Ebondlen + 20.0*old_fene; 
                        // torsion angle
                        float old_tor_fene=0.0;
                        float old_tor_E=0.0;
                        old_tor_E = energyrama(decstrbondlen,pnum,old_tor_fene,rlogduke,ramaduke); 
                        old_tor_E = 4.00*old_tor_E + old_tor_fene;   

                        float E_dist = 0.0;
                        for(int iu=0;iu<N_dist;iu++)
                        {
                            int ik = dist_pre[iu].ires;
                            int jk = dist_pre[iu].jres;
                            float tmp_x = dist_pre[iu].dv;
                            float tmp_y = Distance_point(tmp_mat[3*ik+1].x_,tmp_mat[3*jk+1].x_);
                    //        float tmp_z=sqrt((tmp_x-tmp_y)*(tmp_x-tmp_y))-dist_pre[iu].ave;
                            float tmp_z=tmp_y-dist_pre[iu].ave;
                            E_dist = E_dist + (log(tmp_z*tmp_z+1))/sqrt(dist_pre[iu].stdx);

                        }         
                        E_dist=E_dist/float(pnum*(rand_b2-rand_a2+1)); 

                        float E_con = 0.0;
                        for(int iu=0;iu<N_distg;iu++)
                        {
                            int ik = dist_map[iu].ires;
                            int jk = dist_map[iu].jres;
                            float tmp_x = dist_map[iu].dv;
                            float tmp_y = Distance_point(tmp_mat[3*ik+1].x_,tmp_mat[3*jk+1].x_);
                            E_con = E_con + sqrt((tmp_x-tmp_y)*(tmp_x-tmp_y));
                        }            
                        E_con=E_con/float((rand_b2-rand_a2+1));  

                        float segdisteng0=0.0;
                        if(TF_segmx)
                        {
                            int segint=segm.size();
                            int ttu=1;
                            for(int il=0;il<segint;il++)
                            {
                        //        int intu=segnum[il];
                                int il0=segm[il].init;
                                int iln=segm[il].term;
                                for(int iul=il0;iul<iln;iul++)
                                {
                                    segdisteng0=segdisteng0+Distance_point(tmp_mat[3*iul+1].x_,pointsBy[3*iul+1].x_);
                                    ttu=ttu+1;
                                }
                            }
                            segdisteng0=segdisteng0/float(ttu);
                            segdisteng0=-segdisteng0;
                    //        segdisteng0=1.0/(segdisteng0+1.0);
                        } 
            /*            if(TF_segmx)
                        {
                            int segint=segnum.size();
                            int ttu=1;
                            for(int il=0;il<segint;il++)
                            {
                                int intu=segnum[il];
                                int il0=segm[intu].init;
                                int iln=segm[intu].term;
                                for(int iul=il0;iul<iln;iul++)
                                {
                                    segdisteng0=segdisteng0+Distance_point(tmp_mat[3*iul+1].x_,pointsBy[3*iul+1].x_);
                                    ttu=ttu+1;
                                }
                            }
                            segdisteng0=segdisteng0/float(ttu);
                            segdisteng0=-segdisteng0;
                    //        segdisteng0=1.0/(segdisteng0+1.0);
                        }            */                                 
                        old_dE = w1*old_CC  + 1.0*old_clash + 1.0*(old_Ehbond)+ old_Ebondlen+ 1.0*old_tor_E + 1.0*E_con + w6*E_dist+segdisteng0+Eres0;// + 0.2*Eres0+ old_Ebondang/500.0;
                   //     new_dE = 100.0*new_CC + 100.0*new_clash/abs(rand_b2-rand_a2+1) + 100.0*(new_Ehbond/abs(pnum))/(1.0-(new_Ehbond/abs(pnum)));
                   //     old_dE = old_CC + old_clash + old_Ehbond/abs(pnum)/(1.0-(old_Ehbond/abs(pnum)));
                   //     cout<< "old_dE: "<<old_dE<<" ";
                   //     cout<<"old size: "<<tmp_mat.size()<<endl;
                   //     cout<<"old coor:  "<<tmp_mat[5].x_[0]<<" "<<tmp_mat[5].x_[1]<<" "<<tmp_mat[5].x_[2]<<endl;

                   //     cout<<"fff"<<endl;

            //            rhoC01 = rhoC0;
                         
                        
                        // Rotate angle
                //        ang=180.0;      
        //                vector<vector<float> > tmm;
        //                for(int tc=0;tc<tmp_mat.size();tc++)
        //                {
        //                    tmm.push_back(tmp_mat[tc].x_);
        //                }
        //                GroupRotation(tmm[0],tmm[tmm.size()-1],angle_rotate,tmm);
        //                for(int tc=0;tc<tmp_mat.size();tc++)
        //                {
        //                    tmp_mat[tc].x_ = tmm[tc];
        //                }                
                //        GroupRotationp(tmp_mat[0].x_,tmp_mat[tmp_mat.size()-1].x_,angle_rotate,tmp_mat);
                        if(randtx>=(dm_e-6)||randtx<=(dm_s+6))
                        {  // fragment in the N terminal or C terminal
                    /*        if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E')
                            {
                                int jmx=0,jnx=0;
                                if(randtx>(pnum-4))
                                {
                                    for(int j=0;j<pnum;j++)
                                    {
                                        jmx=randtx-j;
                                        if(jmx<=0)
                                        {
                                            jmx=0;
                                            break;
                                        }
                                        if(numtoss(bb[jmx].sst)=='C') break;

                                    //    if(numtoss(bb[jmx].sst)!='H'&&numtoss(bb[jmx].sst)!='E') break;
                                    }
                                    randtx = jmx;              
                                }
                                if(randtx<3)
                                {
                                    for(int j=0;j<pnum;j++)
                                    {
                                        jnx=randtx+j;
                                        if(jnx>=pnum-1)
                                        {
                                            jnx=pnum-1;
                                            break;
                                        }
                                        if(numtoss(bb[jnx].sst)=='C') break;
                                    //    if(numtoss(bb[jmx].sst)!='H'&&numtoss(bb[jmx].sst)!='E') break;
                                    }
                                    randtx = jnx;            
                                }                   
                            }  */
        //                    cout<<"kkk"<<endl;
                            // Rotation axis 
                            int rand_movement = randIntCustom(0,1000000)%5;
                        //    rand_movement = 1;
                            if(rand_movement==0)
                            {
                                float asin_theta=2.0*randf0and1()-1.0;
                                float acos_theta=sqrt(1.0-asin_theta*asin_theta);
                                float apha=2.0*PI*randf0and1();
                                float awx=acos_theta*cos(apha);
                                float awy=acos_theta*sin(apha);
                                float awz=asin_theta;
                                // Translation Vector
                                float t0=0.3;
                                float t1=(randf0and1()*2.0-1.0)*t0+0.0;
                                float t2=(randf0and1()*2.0-1.0)*t0+0.0;
                                float t3=(randf0and1()*2.0-1.0)*t0+0.0;

                                // Rotation matrix
                                float anggg=90.0;
                                float angle_rotategg=(2.0*randf0and1()-1.0)*anggg; // rotate angle
                                float asin=sin(angle_rotategg*(PI/180.0));
                                float acos=cos(angle_rotategg*(PI/180.0));
                                float u[3][3];
                                u[0][0]=acos+awx*awx*(1.0-acos);
                                u[0][1]=awx*awy*(1.0-acos)-awz*asin;
                                u[0][2]=awx*awz*(1.0-acos)+awy*asin;
                                u[1][0]=awx*awy*(1.0-acos)+awz*asin;
                                u[1][1]=acos+awy*awy*(1.0-acos);
                                u[1][2]=awy*awz*(1.0-acos)-awx*asin;
                                u[2][0]=awx*awz*(1.0-acos)-awy*asin;
                                u[2][1]=awy*awz*(1.0-acos)+awx*asin;
                                u[2][2]=acos+awz*awz*(1.0-acos);            
                                // Rotation points
                        //      int rand_point = rand()%(pnumx);
                        //      int rand_point = rand()%(pnum);
                        //      axyz[0]=pointsx[rand_point].ca[0];
                        //      axyz[1]=pointsx[rand_point].ca[1];
                        //      axyz[2]=pointsx[rand_point].ca[2];
                                int rand_point = rand_a;
                        //        axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
                        //        axyz[1]=pointsBx[3*rand_point+1].x_[1];
                        //        axyz[2]=pointsBx[3*rand_point+1].x_[2];     

                                int tmp3=0;
                                vector<poseCoord > fin_mat;
                                vector<poseCoord>().swap(fin_mat);
                                fin_mat = tmp_mat;
                                if(randtx>=(dm_e-6))
                                {
                                    rand_point = rand_a;
                                    axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
                                    axyz[1]=pointsBx[3*rand_point+1].x_[1];
                                    axyz[2]=pointsBx[3*rand_point+1].x_[2];

                                    for(int j=3*rand_a+1;j<=3*rand_b+2;j++) // not first N
                                    {
                                        fin_mat[j].x_[0]=t1+axyz[0]+(tmp_mat[j].x_[0]-axyz[0])*u[0][0]+(tmp_mat[j].x_[1]-axyz[1])*u[0][1]+(tmp_mat[j].x_[2]-axyz[2])*u[0][2];
                                        fin_mat[j].x_[1]=t2+axyz[1]+(tmp_mat[j].x_[0]-axyz[0])*u[1][0]+(tmp_mat[j].x_[1]-axyz[1])*u[1][1]+(tmp_mat[j].x_[2]-axyz[2])*u[1][2];
                                        fin_mat[j].x_[2]=t3+axyz[2]+(tmp_mat[j].x_[0]-axyz[0])*u[2][0]+(tmp_mat[j].x_[1]-axyz[1])*u[2][1]+(tmp_mat[j].x_[2]-axyz[2])*u[2][2];
                                        tmp3=tmp3+1;
                                    }                                                  
                                } else
                                {
                                    rand_point = rand_b;
                                    axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
                                    axyz[1]=pointsBx[3*rand_point+1].x_[1];
                                    axyz[2]=pointsBx[3*rand_point+1].x_[2];

                                    for(int j=3*rand_a;j<=3*rand_b+1;j++) // not last N
                                    {
                                        fin_mat[j].x_[0]=t1+axyz[0]+(tmp_mat[j].x_[0]-axyz[0])*u[0][0]+(tmp_mat[j].x_[1]-axyz[1])*u[0][1]+(tmp_mat[j].x_[2]-axyz[2])*u[0][2];
                                        fin_mat[j].x_[1]=t2+axyz[1]+(tmp_mat[j].x_[0]-axyz[0])*u[1][0]+(tmp_mat[j].x_[1]-axyz[1])*u[1][1]+(tmp_mat[j].x_[2]-axyz[2])*u[1][2];
                                        fin_mat[j].x_[2]=t3+axyz[2]+(tmp_mat[j].x_[0]-axyz[0])*u[2][0]+(tmp_mat[j].x_[1]-axyz[1])*u[2][1]+(tmp_mat[j].x_[2]-axyz[2])*u[2][2];
                                        tmp3=tmp3+1;
                                    }                        
                                }
                                vector<poseCoord>().swap(tmp_mat);
                                tmp_mat=fin_mat;
                            } 
                            if(rand_movement==1)
                            {
                                // Rotation axis 
                                float asin_theta=2.0*randf0and1()-1.0;
                                float acos_theta=sqrt(1.0-asin_theta*asin_theta);
                                float apha=2.0*PI*randf0and1();
                                float awx=acos_theta*cos(apha);
                                float awy=acos_theta*sin(apha);
                                float awz=asin_theta;
                                // Translation Vector
                                float t0=0.3;
                                float t1=(randf0and1()*2.0-1.0)*t0+0.0;
                                float t2=(randf0and1()*2.0-1.0)*t0+0.0;
                                float t3=(randf0and1()*2.0-1.0)*t0+0.0;

                                // Rotation matrix
                                float anggg=90.0;
                                float angle_rotategg=(2.0*randf0and1()-1.0)*anggg; // rotate angle
                                float asin=sin(angle_rotategg*(PI/180.0));
                                float acos=cos(angle_rotategg*(PI/180.0));
                                float u[3][3];
                                u[0][0]=acos+awx*awx*(1.0-acos);
                                u[0][1]=awx*awy*(1.0-acos)-awz*asin;
                                u[0][2]=awx*awz*(1.0-acos)+awy*asin;
                                u[1][0]=awx*awy*(1.0-acos)+awz*asin;
                                u[1][1]=acos+awy*awy*(1.0-acos);
                                u[1][2]=awy*awz*(1.0-acos)-awx*asin;
                                u[2][0]=awx*awz*(1.0-acos)-awy*asin;
                                u[2][1]=awy*awz*(1.0-acos)+awx*asin;
                                u[2][2]=acos+awz*awz*(1.0-acos);            
                                // Rotation points
                        //      int rand_point = rand()%(pnumx);
                        //      int rand_point = rand()%(pnum);
                        //      axyz[0]=pointsx[rand_point].ca[0];
                        //      axyz[1]=pointsx[rand_point].ca[1];
                        //      axyz[2]=pointsx[rand_point].ca[2];
                                vector<float> coor_pdb(3,0.0);   
                                for(int j=rand_a;j<=rand_b;j++)
                                {
                                    coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
                                    coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
                                    coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];             
                                }
                                coor_pdb[0] = coor_pdb[0]/abs(rand_b-rand_a+1);
                                coor_pdb[1] = coor_pdb[1]/abs(rand_b-rand_a+1);
                                coor_pdb[2] = coor_pdb[2]/abs(rand_b-rand_a+1);                                                   
                        //        int rand_point = rand_a;
                                axyz[0]= coor_pdb[0];  
                                axyz[1]= coor_pdb[1];
                                axyz[2]= coor_pdb[2];               
                                int tmp3=0;
                                vector<poseCoord > fin_mat;
                            //    fin_mat = tmp_mat;
                                vector<poseCoord>().swap(fin_mat);
                                fin_mat = tmp_mat;
                        /*        for(int j=rand_a+2;j<=rand_b-2;j++)
                                {
                                    fin_mat.push_back(pointsBx[3*j+0]);
                                    fin_mat.push_back(pointsBx[3*j+1]);
                                    fin_mat.push_back(pointsBx[3*j+2]);               
                                    tmp3=tmp3+1;
                                }        */
                            //    old_dE = old_clash + olddis_p_m;                      
                                vector<poseCoord > fin_matx;
                                vector<poseCoord >().swap(fin_matx);
                                fin_matx=fin_mat;
                                tmp3=0;
                                for(int j=3*rand_a;j<=3*rand_b+2;j++) // not last N
                                {
                                    fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                                    fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                                    fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                                    tmp3=tmp3+1;
                                }

                                int numseq=pnum;
                            //    if((rand_b-rand_a-3)<2) continue;
                            //    if(rand_a==0) continue;
                                vector<point3f> decstr(numseq);
                                for(int j=0;j<pnum;j++)
                                {
                                    decstr[j].ss2=numtoss(bb[j].sst);
                                    decstr[j].ssm=numtoss(bb[j].sst);
                                    decstr[j].x=fin_mat[3*j+1].x_[0];
                                    decstr[j].y=fin_mat[3*j+1].x_[1];
                                    decstr[j].z=fin_mat[3*j+1].x_[2];
                                    decstr[j].ptn.x=fin_mat[3*j+0].x_[0];
                                    decstr[j].ptn.y=fin_mat[3*j+0].x_[1];
                                    decstr[j].ptn.z=fin_mat[3*j+0].x_[2];
                                    decstr[j].ptc.x=fin_mat[3*j+2].x_[0];
                                    decstr[j].ptc.y=fin_mat[3*j+2].x_[1];
                                    decstr[j].ptc.z=fin_mat[3*j+2].x_[2];
                    //                tmp_tt=tmp_tt+1;
                                }                        
                                if(randtx>=(dm_e-6))
                                {
                                    if((rand_a2)<1||(rand_a+2)>(pnum-1)) continue;
                                    mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+2);
                                } else
                                {
                                    if((rand_b-2)<1||(rand_b2)>(pnum-1)) continue;
                                    mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b2);
                                }                                                                         
                                tmp_tt=0;
                                int tmp_ttx=0;
                                for(int j=0;j<pnum;j++)
                                {
                                    tmp_mat[3*j+1].x_[0]=decstr[j].x;
                                    tmp_mat[3*j+1].x_[1]=decstr[j].y;
                                    tmp_mat[3*j+1].x_[2]=decstr[j].z;
                                    tmp_mat[3*j+0].x_[0]=decstr[j].ptn.x;
                                    tmp_mat[3*j+0].x_[1]=decstr[j].ptn.y;
                                    tmp_mat[3*j+0].x_[2]=decstr[j].ptn.z;
                                    tmp_mat[3*j+2].x_[0]=decstr[j].ptc.x;
                                    tmp_mat[3*j+2].x_[1]=decstr[j].ptc.y;
                                    tmp_mat[3*j+2].x_[2]=decstr[j].ptc.z;
                                    tmp_tt=tmp_tt+1;
                                }     
                                vector<point3f>().swap(decstr);                            
                            }
                            if(rand_movement==2)
                            {
                        /*        bool flagok0=false;
                                for(int j=rand_a;j<=rand_b;j++)
                                {
                                    if(numtoss(bb[j].sst)!='C') flagok0 = true;
                                }
                                if(flagok0) continue; */
                    //            if((rand_b+1)>(pnum-1)) continue;
                                vector<float> sftv(3,0.0);
                                sftv[0] = (2.0*randf0and1()-1.0);
                                sftv[1] = (2.0*randf0and1()-1.0);
                                sftv[2] = (2.0*randf0and1()-1.0);
                                for(int j=rand_a;j<=rand_b;j++)
                                {
                                    tmp_mat[3*j+0].x_[0] = tmp_mat[3*(j)+0].x_[0] + sftv[0];
                                    tmp_mat[3*j+0].x_[1] = tmp_mat[3*(j)+0].x_[1] + sftv[1];
                                    tmp_mat[3*j+0].x_[2] = tmp_mat[3*(j)+0].x_[2] + sftv[2];
                                    tmp_mat[3*j+1].x_[0] = tmp_mat[3*(j)+1].x_[0] + sftv[0];
                                    tmp_mat[3*j+1].x_[1] = tmp_mat[3*(j)+1].x_[1] + sftv[1];
                                    tmp_mat[3*j+1].x_[2] = tmp_mat[3*(j)+1].x_[2] + sftv[2];
                                    tmp_mat[3*j+2].x_[0] = tmp_mat[3*(j)+2].x_[0] + sftv[0];
                                    tmp_mat[3*j+2].x_[1] = tmp_mat[3*(j)+2].x_[1] + sftv[1];
                                    tmp_mat[3*j+2].x_[2] = tmp_mat[3*(j)+2].x_[2] + sftv[2];
                                }
                                int numseq=pnum;
                            //    if((rand_b-rand_a-3)<2) continue;
                            //    if(rand_a==0) continue;
                                vector<point3f> decstr(numseq);

                                for(int j=0;j<pnum;j++)
                                {
                                    decstr[j].ss2=numtoss(bb[j].sst);
                                    decstr[j].ssm=numtoss(bb[j].sst);
                                    decstr[j].x=tmp_mat[3*j+1].x_[0];
                                    decstr[j].y=tmp_mat[3*j+1].x_[1];
                                    decstr[j].z=tmp_mat[3*j+1].x_[2];
                                    decstr[j].ptn.x=tmp_mat[3*j+0].x_[0];
                                    decstr[j].ptn.y=tmp_mat[3*j+0].x_[1];
                                    decstr[j].ptn.z=tmp_mat[3*j+0].x_[2];
                                    decstr[j].ptc.x=tmp_mat[3*j+2].x_[0];
                                    decstr[j].ptc.y=tmp_mat[3*j+2].x_[1];
                                    decstr[j].ptc.z=tmp_mat[3*j+2].x_[2];
                    //                tmp_tt=tmp_tt+1;
                                }

                                if(randtx>=(dm_e-6))
                                {
                                    if((rand_a2)<1||(rand_a+2)>(pnum-1)) continue;
                                    mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+2);
                                }else
                                {
                                    if((rand_b2)>(pnum-1)||(rand_b-2)<1) continue;
                                    mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b2); 
                                }
                                tmp_tt=0;
                                int tmp_ttx=0;
                                for(int j=0;j<pnum;j++)
                                {
                                    tmp_mat[3*j+1].x_[0]=decstr[j].x;
                                    tmp_mat[3*j+1].x_[1]=decstr[j].y;
                                    tmp_mat[3*j+1].x_[2]=decstr[j].z;
                                    tmp_mat[3*j+0].x_[0]=decstr[j].ptn.x;
                                    tmp_mat[3*j+0].x_[1]=decstr[j].ptn.y;
                                    tmp_mat[3*j+0].x_[2]=decstr[j].ptn.z;
                                    tmp_mat[3*j+2].x_[0]=decstr[j].ptc.x;
                                    tmp_mat[3*j+2].x_[1]=decstr[j].ptc.y;
                                    tmp_mat[3*j+2].x_[2]=decstr[j].ptc.z;
                                    tmp_tt=tmp_tt+1;
                                }     
                                vector<point3f>().swap(decstr);                
                            } 
                            if(rand_movement==6)
                            {
                                int numseq=pnum;
                            //    if((rand_b-rand_a-3)<=2) continue;
                                vector<point3f> decstr(numseq);

                                int tmp_tt=0;
                                int tmp_ttx=0;
                                for(int j=0;j<pnum;j++)
                                {
                                    decstr[j].ss2=numtoss(bb[j].sst);
                                    decstr[j].ssm=numtoss(bb[j].sst);
                                    decstr[j].x=tmp_mat[3*tmp_tt+1].x_[0];
                                    decstr[j].y=tmp_mat[3*tmp_tt+1].x_[1];
                                    decstr[j].z=tmp_mat[3*tmp_tt+1].x_[2];
                                    decstr[j].ptn.x=tmp_mat[3*tmp_tt+0].x_[0];
                                    decstr[j].ptn.y=tmp_mat[3*tmp_tt+0].x_[1];
                                    decstr[j].ptn.z=tmp_mat[3*tmp_tt+0].x_[2];
                                    decstr[j].ptc.x=tmp_mat[3*tmp_tt+2].x_[0];
                                    decstr[j].ptc.y=tmp_mat[3*tmp_tt+2].x_[1];
                                    decstr[j].ptc.z=tmp_mat[3*tmp_tt+2].x_[2];
                                    tmp_tt=tmp_tt+1;
                                }
                                if(rand_a==0) continue;
                                mcfragsweepLMP2(decstr,numseq,rand_a,rand_b);
                                tmp_tt=0;
                                tmp_ttx=0;
                                for(int j=0;j<pnum;j++)
                                {
                                    tmp_mat[3*tmp_tt+1].x_[0]=decstr[j].x;
                                    tmp_mat[3*tmp_tt+1].x_[1]=decstr[j].y;
                                    tmp_mat[3*tmp_tt+1].x_[2]=decstr[j].z;
                                    tmp_mat[3*tmp_tt+0].x_[0]=decstr[j].ptn.x;
                                    tmp_mat[3*tmp_tt+0].x_[1]=decstr[j].ptn.y;
                                    tmp_mat[3*tmp_tt+0].x_[2]=decstr[j].ptn.z;
                                    tmp_mat[3*tmp_tt+2].x_[0]=decstr[j].ptc.x;
                                    tmp_mat[3*tmp_tt+2].x_[1]=decstr[j].ptc.y;
                                    tmp_mat[3*tmp_tt+2].x_[2]=decstr[j].ptc.z;
                                    tmp_tt=tmp_tt+1;
                                }
                                vector<point3f>().swap(decstr);                
                            } 
                            if(rand_movement==3) //rotate two teminal around medium point
                            {
                                int randd = int((rand_a + rand_b)/2);
                                float asin_theta=2.0*randf0and1()-1.0;
                                float acos_theta=sqrt(1.0-asin_theta*asin_theta);
                                float apha=2.0*PI*randf0and1();
                                float awx=acos_theta*cos(apha);
                                float awy=acos_theta*sin(apha);
                                float awz=asin_theta;
                                // Translation Vector
                                float t0=0.3;
                                float t1=(randf0and1()*2.0-1.0)*t0+0.0;
                                float t2=(randf0and1()*2.0-1.0)*t0+0.0;
                                float t3=(randf0and1()*2.0-1.0)*t0+0.0;

                                // Rotation matrix
                                float anggg=90.0;
                                float angle_rotategg=(2.0*randf0and1()-1.0)*anggg; // rotate angle
                                float asin=sin(angle_rotategg*(PI/180.0));
                                float acos=cos(angle_rotategg*(PI/180.0));
                                float u[3][3];
                                u[0][0]=acos+awx*awx*(1.0-acos);
                                u[0][1]=awx*awy*(1.0-acos)-awz*asin;
                                u[0][2]=awx*awz*(1.0-acos)+awy*asin;
                                u[1][0]=awx*awy*(1.0-acos)+awz*asin;
                                u[1][1]=acos+awy*awy*(1.0-acos);
                                u[1][2]=awy*awz*(1.0-acos)-awx*asin;
                                u[2][0]=awx*awz*(1.0-acos)-awy*asin;
                                u[2][1]=awy*awz*(1.0-acos)+awx*asin;
                                u[2][2]=acos+awz*awz*(1.0-acos); 

                                axyz[0]=pointsBx[3*randd+1].x_[0]; // CA atom
                                axyz[1]=pointsBx[3*randd+1].x_[1];
                                axyz[2]=pointsBx[3*randd+1].x_[2];
                                int tmp3=0;
                                vector<poseCoord > fin_mat;
                            //    fin_mat = tmp_mat;
                                vector<poseCoord>().swap(fin_mat);
                                fin_mat = tmp_mat;  
                                vector<poseCoord > fin_matx;
                                vector<poseCoord >().swap(fin_matx);
                                fin_matx=fin_mat;
                                tmp3=0;
                                for(int j=3*rand_a;j<=3*randd+1;j++) // not last N
                                {
                                    fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                                    fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                                    fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                                    tmp3=tmp3+1;
                                }
                                tmp3=0;
                                for(int j=3*randd+1;j<=3*rand_b+2;j++) // not last N
                                {
                                    fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                                    fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                                    fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                                    tmp3=tmp3+1;
                                }
                                int numseq=pnum;
                            //    if((rand_b-rand_a-3)<2) continue;
                            //    if(rand_a==0) continue;
                                vector<point3f> decstr(numseq);
                                for(int j=0;j<pnum;j++)
                                {
                                    decstr[j].ss2=numtoss(bb[j].sst);
                                    decstr[j].ssm=numtoss(bb[j].sst);
                                    decstr[j].x=fin_mat[3*j+1].x_[0];
                                    decstr[j].y=fin_mat[3*j+1].x_[1];
                                    decstr[j].z=fin_mat[3*j+1].x_[2];
                                    decstr[j].ptn.x=fin_mat[3*j+0].x_[0];
                                    decstr[j].ptn.y=fin_mat[3*j+0].x_[1];
                                    decstr[j].ptn.z=fin_mat[3*j+0].x_[2];
                                    decstr[j].ptc.x=fin_mat[3*j+2].x_[0];
                                    decstr[j].ptc.y=fin_mat[3*j+2].x_[1];
                                    decstr[j].ptc.z=fin_mat[3*j+2].x_[2];
                    //                tmp_tt=tmp_tt+1;
                                }
                                if(randtx>=(dm_e-6))
                                {
                                    if((rand_a2)<1||(rand_a+2)>(pnum-1)) continue;
                                    mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+2); 
                                } else
                                {
                                    if((rand_b2)>(pnum-1)||(rand_b-2)<1) continue;
                                    mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b2);
                                }
                     
                                tmp_tt=0;
                                int tmp_ttx=0;
                                for(int j=0;j<pnum;j++)
                                {
                                    tmp_mat[3*j+1].x_[0]=decstr[j].x;
                                    tmp_mat[3*j+1].x_[1]=decstr[j].y;
                                    tmp_mat[3*j+1].x_[2]=decstr[j].z;
                                    tmp_mat[3*j+0].x_[0]=decstr[j].ptn.x;
                                    tmp_mat[3*j+0].x_[1]=decstr[j].ptn.y;
                                    tmp_mat[3*j+0].x_[2]=decstr[j].ptn.z;
                                    tmp_mat[3*j+2].x_[0]=decstr[j].ptc.x;
                                    tmp_mat[3*j+2].x_[1]=decstr[j].ptc.y;
                                    tmp_mat[3*j+2].x_[2]=decstr[j].ptc.z;
                                    tmp_tt=tmp_tt+1;
                                }     
                                vector<point3f>().swap(decstr);                                                                   
                            }
                            if(rand_movement==4)
                            {
                                angle_rotate = (2.0*randf0and1()-1.0)*180.0; // rotate angle 
                                GroupRotationpid(tmp_mat[3*rand_a+1].x_,tmp_mat[3*rand_b+2-1].x_,angle_rotate,tmp_mat,3*rand_a+1,3*rand_b+1);
                            }                                                                                                                                    
                           
                        } else 
                        {
                            int rand_movement = randIntCustom(0,1000000)%7;
                    //        rand_movement = 7;
        //                    cout<<"sssssssssssssssss: "<<rand_movement<<endl;
                            if(rand_movement==0)
                            {
                                angle_rotate = (2.0*randf0and1()-1.0)*90.0; // rotate angle 
                                GroupRotationpid(tmp_mat[3*rand_a+1].x_,tmp_mat[3*rand_b+2-1].x_,angle_rotate,tmp_mat,3*rand_a+1,3*rand_b+1);
                            }                  
                            if(rand_movement==1)
                            {
                        /*        bool flagok0=false;
                                for(int j=rand_a;j<=rand_b;j++)
                                {
                                    if(numtoss(bb[j].sst)!='C') flagok0 = true;
                                }
                                if(flagok0) continue; */
                    //            if((rand_b+1)>(pnum-1)) continue;
                                vector<float> sftv(3,0.0);
                                sftv[0] = (2.0*randf0and1()-1.0);
                                sftv[1] = (2.0*randf0and1()-1.0);
                                sftv[2] = (2.0*randf0and1()-1.0);
                                for(int j=rand_a;j<=rand_b;j++)
                                {
                                    tmp_mat[3*j+0].x_[0] = tmp_mat[3*(j)+0].x_[0] + sftv[0];
                                    tmp_mat[3*j+0].x_[1] = tmp_mat[3*(j)+0].x_[1] + sftv[1];
                                    tmp_mat[3*j+0].x_[2] = tmp_mat[3*(j)+0].x_[2] + sftv[2];
                                    tmp_mat[3*j+1].x_[0] = tmp_mat[3*(j)+1].x_[0] + sftv[0];
                                    tmp_mat[3*j+1].x_[1] = tmp_mat[3*(j)+1].x_[1] + sftv[1];
                                    tmp_mat[3*j+1].x_[2] = tmp_mat[3*(j)+1].x_[2] + sftv[2];
                                    tmp_mat[3*j+2].x_[0] = tmp_mat[3*(j)+2].x_[0] + sftv[0];
                                    tmp_mat[3*j+2].x_[1] = tmp_mat[3*(j)+2].x_[1] + sftv[1];
                                    tmp_mat[3*j+2].x_[2] = tmp_mat[3*(j)+2].x_[2] + sftv[2];
                                }
                                int numseq=pnum;
                            //    if((rand_b-rand_a-3)<2) continue;
                            //    if(rand_a==0) continue;
                                vector<point3f> decstr(numseq);

                                for(int j=0;j<pnum;j++)
                                {
                                    decstr[j].ss2=numtoss(bb[j].sst);
                                    decstr[j].ssm=numtoss(bb[j].sst);
                                    decstr[j].x=tmp_mat[3*j+1].x_[0];
                                    decstr[j].y=tmp_mat[3*j+1].x_[1];
                                    decstr[j].z=tmp_mat[3*j+1].x_[2];
                                    decstr[j].ptn.x=tmp_mat[3*j+0].x_[0];
                                    decstr[j].ptn.y=tmp_mat[3*j+0].x_[1];
                                    decstr[j].ptn.z=tmp_mat[3*j+0].x_[2];
                                    decstr[j].ptc.x=tmp_mat[3*j+2].x_[0];
                                    decstr[j].ptc.y=tmp_mat[3*j+2].x_[1];
                                    decstr[j].ptc.z=tmp_mat[3*j+2].x_[2];
                    //                tmp_tt=tmp_tt+1;
                                }
                                if((rand_a-2)<1||(rand_b-2)<1) continue;
                                if((rand_b+2)>(pnum-1)||(rand_a+2)>(pnum-1)) continue;
                                mcfragsweepLMP2(decstr,numseq,rand_a-2,rand_a+2);
                                mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b+2); 
                                tmp_tt=0;
                                int tmp_ttx=0;
                                for(int j=0;j<pnum;j++)
                                {
                                    tmp_mat[3*j+1].x_[0]=decstr[j].x;
                                    tmp_mat[3*j+1].x_[1]=decstr[j].y;
                                    tmp_mat[3*j+1].x_[2]=decstr[j].z;
                                    tmp_mat[3*j+0].x_[0]=decstr[j].ptn.x;
                                    tmp_mat[3*j+0].x_[1]=decstr[j].ptn.y;
                                    tmp_mat[3*j+0].x_[2]=decstr[j].ptn.z;
                                    tmp_mat[3*j+2].x_[0]=decstr[j].ptc.x;
                                    tmp_mat[3*j+2].x_[1]=decstr[j].ptc.y;
                                    tmp_mat[3*j+2].x_[2]=decstr[j].ptc.z;
                                    tmp_tt=tmp_tt+1;
                                }     
                                vector<point3f>().swap(decstr);                
                            }                     
                            if(rand_movement==2) //wrong
                            {
                                // Rotation axis 
                                float asin_theta=2.0*randf0and1()-1.0;
                                float acos_theta=sqrt(1.0-asin_theta*asin_theta);
                                float apha=2.0*PI*randf0and1();
                                float awx=acos_theta*cos(apha);
                                float awy=acos_theta*sin(apha);
                                float awz=asin_theta;
                                // Translation Vector
                                float t0=0.3;
                                float t1=(randf0and1()*2.0-1.0)*t0+0.0;
                                float t2=(randf0and1()*2.0-1.0)*t0+0.0;
                                float t3=(randf0and1()*2.0-1.0)*t0+0.0;

                                // Rotation matrix
                                float anggg=180.0;
                                float angle_rotategg=(2.0*randf0and1()-1.0)*anggg; // rotate angle
                                float asin=sin(angle_rotategg*(PI/180.0));
                                float acos=cos(angle_rotategg*(PI/180.0));
                                float u[3][3];
                                u[0][0]=acos+awx*awx*(1.0-acos);
                                u[0][1]=awx*awy*(1.0-acos)-awz*asin;
                                u[0][2]=awx*awz*(1.0-acos)+awy*asin;
                                u[1][0]=awx*awy*(1.0-acos)+awz*asin;
                                u[1][1]=acos+awy*awy*(1.0-acos);
                                u[1][2]=awy*awz*(1.0-acos)-awx*asin;
                                u[2][0]=awx*awz*(1.0-acos)-awy*asin;
                                u[2][1]=awy*awz*(1.0-acos)+awx*asin;
                                u[2][2]=acos+awz*awz*(1.0-acos);            
                                // Rotation points
                        //      int rand_point = rand()%(pnumx);
                        //      int rand_point = rand()%(pnum);
                        //      axyz[0]=pointsx[rand_point].ca[0];
                        //      axyz[1]=pointsx[rand_point].ca[1];
                        //      axyz[2]=pointsx[rand_point].ca[2];
                        /*        vector<float> coor_pdb(3,0.0);   
                                for(int j=rand_a;j<=rand_b;j++)
                                {
                                    coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
                                    coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
                                    coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];             
                                }
                                coor_pdb[0] = coor_pdb[0]/abs(rand_b-rand_a+1);
                                coor_pdb[1] = coor_pdb[1]/abs(rand_b-rand_a+1);
                                coor_pdb[2] = coor_pdb[2]/abs(rand_b-rand_a+1);                                                   
                        //        int rand_point = rand_a;
                                axyz[0]= coor_pdb[0];  
                                axyz[1]= coor_pdb[1];
                                axyz[2]= coor_pdb[2];               */
                                float rand_mov0=(randf0and1()*2.0-1.0);
                                if(rand_mov0>0)
                                {
                                    axyz[0]=pointsBx[3*rand_a+1].x_[0]; // CA atom
                                    axyz[1]=pointsBx[3*rand_a+1].x_[1];
                                    axyz[2]=pointsBx[3*rand_a+1].x_[2];
                                }
                                else 
                                {
                                    axyz[0]=pointsBx[3*rand_b+1].x_[0]; // CA atom
                                    axyz[1]=pointsBx[3*rand_b+1].x_[1];
                                    axyz[2]=pointsBx[3*rand_b+1].x_[2];                            
                                }
                                int tmp3=0;
                                vector<poseCoord > fin_mat;
                            //    fin_mat = tmp_mat;
                                vector<poseCoord>().swap(fin_mat);
                                fin_mat = tmp_mat;
                        /*        for(int j=rand_a+2;j<=rand_b-2;j++)
                                {
                                    fin_mat.push_back(pointsBx[3*j+0]);
                                    fin_mat.push_back(pointsBx[3*j+1]);
                                    fin_mat.push_back(pointsBx[3*j+2]);               
                                    tmp3=tmp3+1;
                                }        */
                            //    old_dE = old_clash + olddis_p_m;                      
                                vector<poseCoord > fin_matx;
                                vector<poseCoord >().swap(fin_matx);
                                fin_matx=fin_mat;
                                tmp3=0;
                                if(rand_mov0>0)
                                {
                                    for(int j=3*rand_a+1;j<=3*rand_b+2;j++) // not last N
                                    {
                                        fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                                        fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                                        fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                                        tmp3=tmp3+1;
                                    }
                                }else
                                {
                                    for(int j=3*rand_a;j<=3*rand_b+1;j++) // not last N
                                    {
                                        fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                                        fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                                        fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                                        tmp3=tmp3+1;
                                    }
                                }

                                int numseq=pnum;
                            //    if((rand_b-rand_a-3)<2) continue;
                            //    if(rand_a==0) continue;
                                vector<point3f> decstr(numseq);
                                for(int j=0;j<pnum;j++)
                                {
                                    decstr[j].ss2=numtoss(bb[j].sst);
                                    decstr[j].ssm=numtoss(bb[j].sst);
                                    decstr[j].x=fin_mat[3*j+1].x_[0];
                                    decstr[j].y=fin_mat[3*j+1].x_[1];
                                    decstr[j].z=fin_mat[3*j+1].x_[2];
                                    decstr[j].ptn.x=fin_mat[3*j+0].x_[0];
                                    decstr[j].ptn.y=fin_mat[3*j+0].x_[1];
                                    decstr[j].ptn.z=fin_mat[3*j+0].x_[2];
                                    decstr[j].ptc.x=fin_mat[3*j+2].x_[0];
                                    decstr[j].ptc.y=fin_mat[3*j+2].x_[1];
                                    decstr[j].ptc.z=fin_mat[3*j+2].x_[2];
                    //                tmp_tt=tmp_tt+1;
                                }
                                if((rand_a-2)<1||(rand_b-2)<1) continue;
                                if((rand_b+2)>(pnum-1)||(rand_a+2)>(pnum-1)) continue;
                                if(rand_mov0>0)
                                {
                                    vector<int> beg_end(2,0);
                                    getLinLen(rand_b+1,fin_mat,bb,beg_end);
                                    mcfragsweepLMP2(decstr,numseq,beg_end[0],beg_end[1]);
                                } else
                                {
                                    vector<int> beg_end(2,0);
                                    getLinLen(rand_a,fin_mat,bb,beg_end);
                                    mcfragsweepLMP2(decstr,numseq,beg_end[0],beg_end[1]);
                                }
                                tmp_tt=0;
                                int tmp_ttx=0;
                                for(int j=0;j<pnum;j++)
                                {
                                    tmp_mat[3*j+1].x_[0]=decstr[j].x;
                                    tmp_mat[3*j+1].x_[1]=decstr[j].y;
                                    tmp_mat[3*j+1].x_[2]=decstr[j].z;
                                    tmp_mat[3*j+0].x_[0]=decstr[j].ptn.x;
                                    tmp_mat[3*j+0].x_[1]=decstr[j].ptn.y;
                                    tmp_mat[3*j+0].x_[2]=decstr[j].ptn.z;
                                    tmp_mat[3*j+2].x_[0]=decstr[j].ptc.x;
                                    tmp_mat[3*j+2].x_[1]=decstr[j].ptc.y;
                                    tmp_mat[3*j+2].x_[2]=decstr[j].ptc.z;
                                    tmp_tt=tmp_tt+1;
                                }     
                                vector<point3f>().swap(decstr);                       
                            } 
                            if(rand_movement==3) //wrong
                            {
                                // Rotation axis 
                                float asin_theta=2.0*randf0and1()-1.0;
                                float acos_theta=sqrt(1.0-asin_theta*asin_theta);
                                float apha=2.0*PI*randf0and1();
                                float awx=acos_theta*cos(apha);
                                float awy=acos_theta*sin(apha);
                                float awz=asin_theta;
                                // Translation Vector
                                float t0=0.3;
                                float t1=(randf0and1()*2.0-1.0)*t0+0.0;
                                float t2=(randf0and1()*2.0-1.0)*t0+0.0;
                                float t3=(randf0and1()*2.0-1.0)*t0+0.0;

                                // Rotation matrix
                                float anggg=90.0;
                                float angle_rotategg=(2.0*randf0and1()-1.0)*anggg; // rotate angle
                                float asin=sin(angle_rotategg*(PI/180.0));
                                float acos=cos(angle_rotategg*(PI/180.0));
                                float u[3][3];
                                u[0][0]=acos+awx*awx*(1.0-acos);
                                u[0][1]=awx*awy*(1.0-acos)-awz*asin;
                                u[0][2]=awx*awz*(1.0-acos)+awy*asin;
                                u[1][0]=awx*awy*(1.0-acos)+awz*asin;
                                u[1][1]=acos+awy*awy*(1.0-acos);
                                u[1][2]=awy*awz*(1.0-acos)-awx*asin;
                                u[2][0]=awx*awz*(1.0-acos)-awy*asin;
                                u[2][1]=awy*awz*(1.0-acos)+awx*asin;
                                u[2][2]=acos+awz*awz*(1.0-acos);            
                                // Rotation points
                        //      int rand_point = rand()%(pnumx);
                        //      int rand_point = rand()%(pnum);
                        //      axyz[0]=pointsx[rand_point].ca[0];
                        //      axyz[1]=pointsx[rand_point].ca[1];
                        //      axyz[2]=pointsx[rand_point].ca[2];
                                vector<float> coor_pdb(3,0.0);   
                                for(int j=rand_a;j<=rand_b;j++)
                                {
                                    coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
                                    coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
                                    coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];             
                                }
                                coor_pdb[0] = coor_pdb[0]/abs(rand_b-rand_a+1);
                                coor_pdb[1] = coor_pdb[1]/abs(rand_b-rand_a+1);
                                coor_pdb[2] = coor_pdb[2]/abs(rand_b-rand_a+1);                                                   
                        //        int rand_point = rand_a;
                                axyz[0]= coor_pdb[0];  
                                axyz[1]= coor_pdb[1];
                                axyz[2]= coor_pdb[2];               
                                int tmp3=0;
                                vector<poseCoord > fin_mat;
                            //    fin_mat = tmp_mat;
                                vector<poseCoord>().swap(fin_mat);
                                fin_mat = tmp_mat;
                        /*        for(int j=rand_a+2;j<=rand_b-2;j++)
                                {
                                    fin_mat.push_back(pointsBx[3*j+0]);
                                    fin_mat.push_back(pointsBx[3*j+1]);
                                    fin_mat.push_back(pointsBx[3*j+2]);               
                                    tmp3=tmp3+1;
                                }        */
                            //    old_dE = old_clash + olddis_p_m;                      
                                vector<poseCoord > fin_matx;
                                vector<poseCoord >().swap(fin_matx);
                                fin_matx=fin_mat;
                                tmp3=0;
                                for(int j=3*rand_a;j<=3*rand_b+2;j++) // not last N
                                {
                                    fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                                    fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                                    fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                                    tmp3=tmp3+1;
                                }

                                int numseq=pnum;
                            //    if((rand_b-rand_a-3)<2) continue;
                            //    if(rand_a==0) continue;
                                vector<point3f> decstr(numseq);
                                for(int j=0;j<pnum;j++)
                                {
                                    decstr[j].ss2=numtoss(bb[j].sst);
                                    decstr[j].ssm=numtoss(bb[j].sst);
                                    decstr[j].x=fin_mat[3*j+1].x_[0];
                                    decstr[j].y=fin_mat[3*j+1].x_[1];
                                    decstr[j].z=fin_mat[3*j+1].x_[2];
                                    decstr[j].ptn.x=fin_mat[3*j+0].x_[0];
                                    decstr[j].ptn.y=fin_mat[3*j+0].x_[1];
                                    decstr[j].ptn.z=fin_mat[3*j+0].x_[2];
                                    decstr[j].ptc.x=fin_mat[3*j+2].x_[0];
                                    decstr[j].ptc.y=fin_mat[3*j+2].x_[1];
                                    decstr[j].ptc.z=fin_mat[3*j+2].x_[2];
                    //                tmp_tt=tmp_tt+1;
                                }
                                if((rand_a-2)<1||(rand_b-2)<1) continue;
                                if((rand_b+2)>(pnum-1)||(rand_a+2)>(pnum-1)) continue;
                                vector<int> beg_end(2,0);
                                getLinLen(rand_b+1,fin_mat,bb,beg_end);
                                mcfragsweepLMP2(decstr,numseq,beg_end[0],beg_end[1]);
                                getLinLen(rand_a,fin_mat,bb,beg_end);                            
                                mcfragsweepLMP2(decstr,numseq,beg_end[0],beg_end[1]); 
                            //    mcfragsweepLMP2(decstr,numseq,rand_a-2,rand_b+2);                    
                                tmp_tt=0;
                                int tmp_ttx=0;
                                for(int j=0;j<pnum;j++)
                                {
                                    tmp_mat[3*j+1].x_[0]=decstr[j].x;
                                    tmp_mat[3*j+1].x_[1]=decstr[j].y;
                                    tmp_mat[3*j+1].x_[2]=decstr[j].z;
                                    tmp_mat[3*j+0].x_[0]=decstr[j].ptn.x;
                                    tmp_mat[3*j+0].x_[1]=decstr[j].ptn.y;
                                    tmp_mat[3*j+0].x_[2]=decstr[j].ptn.z;
                                    tmp_mat[3*j+2].x_[0]=decstr[j].ptc.x;
                                    tmp_mat[3*j+2].x_[1]=decstr[j].ptc.y;
                                    tmp_mat[3*j+2].x_[2]=decstr[j].ptc.z;
                                    tmp_tt=tmp_tt+1;
                                }     
                                vector<point3f>().swap(decstr);                        
                            }
                            if(rand_movement==4) //rotate two teminal around medium point
                            {
                                int randd = int((rand_a + rand_b)/2);
                                float asin_theta=2.0*randf0and1()-1.0;
                                float acos_theta=sqrt(1.0-asin_theta*asin_theta);
                                float apha=2.0*PI*randf0and1();
                                float awx=acos_theta*cos(apha);
                                float awy=acos_theta*sin(apha);
                                float awz=asin_theta;
                                // Translation Vector
                                float t0=0.3;
                                float t1=(randf0and1()*2.0-1.0)*t0+0.0;
                                float t2=(randf0and1()*2.0-1.0)*t0+0.0;
                                float t3=(randf0and1()*2.0-1.0)*t0+0.0;

                                // Rotation matrix
                                float anggg=90.0;
                                float angle_rotategg=(2.0*randf0and1()-1.0)*anggg; // rotate angle
                                float asin=sin(angle_rotategg*(PI/180.0));
                                float acos=cos(angle_rotategg*(PI/180.0));
                                float u[3][3];
                                u[0][0]=acos+awx*awx*(1.0-acos);
                                u[0][1]=awx*awy*(1.0-acos)-awz*asin;
                                u[0][2]=awx*awz*(1.0-acos)+awy*asin;
                                u[1][0]=awx*awy*(1.0-acos)+awz*asin;
                                u[1][1]=acos+awy*awy*(1.0-acos);
                                u[1][2]=awy*awz*(1.0-acos)-awx*asin;
                                u[2][0]=awx*awz*(1.0-acos)-awy*asin;
                                u[2][1]=awy*awz*(1.0-acos)+awx*asin;
                                u[2][2]=acos+awz*awz*(1.0-acos); 

                                axyz[0]=pointsBx[3*randd+1].x_[0]; // CA atom
                                axyz[1]=pointsBx[3*randd+1].x_[1];
                                axyz[2]=pointsBx[3*randd+1].x_[2];
                                int tmp3=0;
                                vector<poseCoord > fin_mat;
                            //    fin_mat = tmp_mat;
                                vector<poseCoord>().swap(fin_mat);
                                fin_mat = tmp_mat;  
                                vector<poseCoord > fin_matx;
                                vector<poseCoord >().swap(fin_matx);
                                fin_matx=fin_mat;
                                tmp3=0;
                                for(int j=3*rand_a;j<=3*randd+1;j++) // not last N
                                {
                                    fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                                    fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                                    fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                                    tmp3=tmp3+1;
                                }
                                tmp3=0;
                                for(int j=3*randd+1;j<=3*rand_b+2;j++) // not last N
                                {
                                    fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                                    fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                                    fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                                    tmp3=tmp3+1;
                                }
                                int numseq=pnum;
                            //    if((rand_b-rand_a-3)<2) continue;
                            //    if(rand_a==0) continue;
                                vector<point3f> decstr(numseq);
                                for(int j=0;j<pnum;j++)
                                {
                                    decstr[j].ss2=numtoss(bb[j].sst);
                                    decstr[j].ssm=numtoss(bb[j].sst);
                                    decstr[j].x=fin_mat[3*j+1].x_[0];
                                    decstr[j].y=fin_mat[3*j+1].x_[1];
                                    decstr[j].z=fin_mat[3*j+1].x_[2];
                                    decstr[j].ptn.x=fin_mat[3*j+0].x_[0];
                                    decstr[j].ptn.y=fin_mat[3*j+0].x_[1];
                                    decstr[j].ptn.z=fin_mat[3*j+0].x_[2];
                                    decstr[j].ptc.x=fin_mat[3*j+2].x_[0];
                                    decstr[j].ptc.y=fin_mat[3*j+2].x_[1];
                                    decstr[j].ptc.z=fin_mat[3*j+2].x_[2];
                    //                tmp_tt=tmp_tt+1;
                                }
                                if((rand_a-2)<1||(rand_b-2)<1) continue;
                                if((rand_b+2)>(pnum-1)||(rand_a+2)>(pnum-1)) continue;
                                mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b+2);                            
                                mcfragsweepLMP2(decstr,numseq,rand_a-2,rand_a+2);                     
                                tmp_tt=0;
                                int tmp_ttx=0;
                                for(int j=0;j<pnum;j++)
                                {
                                    tmp_mat[3*j+1].x_[0]=decstr[j].x;
                                    tmp_mat[3*j+1].x_[1]=decstr[j].y;
                                    tmp_mat[3*j+1].x_[2]=decstr[j].z;
                                    tmp_mat[3*j+0].x_[0]=decstr[j].ptn.x;
                                    tmp_mat[3*j+0].x_[1]=decstr[j].ptn.y;
                                    tmp_mat[3*j+0].x_[2]=decstr[j].ptn.z;
                                    tmp_mat[3*j+2].x_[0]=decstr[j].ptc.x;
                                    tmp_mat[3*j+2].x_[1]=decstr[j].ptc.y;
                                    tmp_mat[3*j+2].x_[2]=decstr[j].ptc.z;
                                    tmp_tt=tmp_tt+1;
                                }     
                                vector<point3f>().swap(decstr);                                                                   
                            } 
                            if(rand_movement==5)
                            {
                                int numseq=pnum;
                            //    if((rand_b-rand_a-3)<=2) continue;
                                vector<point3f> decstr(numseq);

                                int tmp_tt=0;
                                int tmp_ttx=0;
                                for(int j=0;j<pnum;j++)
                                {
                                    decstr[j].ss2=numtoss(bb[j].sst);
                                    decstr[j].ssm=numtoss(bb[j].sst);
                                    decstr[j].x=tmp_mat[3*tmp_tt+1].x_[0];
                                    decstr[j].y=tmp_mat[3*tmp_tt+1].x_[1];
                                    decstr[j].z=tmp_mat[3*tmp_tt+1].x_[2];
                                    decstr[j].ptn.x=tmp_mat[3*tmp_tt+0].x_[0];
                                    decstr[j].ptn.y=tmp_mat[3*tmp_tt+0].x_[1];
                                    decstr[j].ptn.z=tmp_mat[3*tmp_tt+0].x_[2];
                                    decstr[j].ptc.x=tmp_mat[3*tmp_tt+2].x_[0];
                                    decstr[j].ptc.y=tmp_mat[3*tmp_tt+2].x_[1];
                                    decstr[j].ptc.z=tmp_mat[3*tmp_tt+2].x_[2];
                                    tmp_tt=tmp_tt+1;
                                }
                                if(rand_a==0) continue;
                            //    mcfragsweepLMP2(decstr,numseq,rand_a,rand_b);
                                mcmovementLMP(decstr,numseq,rand_a,rand_b);
                                tmp_tt=0;
                                tmp_ttx=0;
                                for(int j=0;j<pnum;j++)
                                {
                                    tmp_mat[3*tmp_tt+1].x_[0]=decstr[j].x;
                                    tmp_mat[3*tmp_tt+1].x_[1]=decstr[j].y;
                                    tmp_mat[3*tmp_tt+1].x_[2]=decstr[j].z;
                                    tmp_mat[3*tmp_tt+0].x_[0]=decstr[j].ptn.x;
                                    tmp_mat[3*tmp_tt+0].x_[1]=decstr[j].ptn.y;
                                    tmp_mat[3*tmp_tt+0].x_[2]=decstr[j].ptn.z;
                                    tmp_mat[3*tmp_tt+2].x_[0]=decstr[j].ptc.x;
                                    tmp_mat[3*tmp_tt+2].x_[1]=decstr[j].ptc.y;
                                    tmp_mat[3*tmp_tt+2].x_[2]=decstr[j].ptc.z;
                                    tmp_tt=tmp_tt+1;
                                }
                                vector<point3f>().swap(decstr);                
                            }
                            if(rand_movement==6) //shift
                            {
                        /*        bool flagok0=false;
                                for(int j=rand_a;j<=rand_b;j++)
                                {
                                    if(numtoss(bb[j].sst)!='C') flagok0 = true;
                                }  
                                if(flagok0) continue; */ 
                                int randg= 1;//randIntCustom(1,2);
                                float rand_mov0=(randf0and1()*2.0-1.0);
                                if(rand_mov0>=0) // left shift
                                {
                                //    if((rand_b+randg)>(pnum-1)) continue;
                                    if((rand_b+randg)>(pnum-1)) continue;
                                    for(int j=rand_a;j<=rand_b;j++)
                                    {
        //                                vector<float> sftv(3,0.0);
        //                                sftv[0]= tmp_mat[3*j+1].x_[0]-tmp_mat[3*(j+randg)+1].x_[0];
        //                                sftv[1]= tmp_mat[3*j+1].x_[1]-tmp_mat[3*(j+randg)+1].x_[1];
        //                                sftv[2]= tmp_mat[3*j+1].x_[2]-tmp_mat[3*(j+randg)+1].x_[2];                                
                                        tmp_mat[3*j+0].x_[0] = tmp_mat[3*(j+randg)+0].x_[0];// + sftv[0];
                                        tmp_mat[3*j+0].x_[1] = tmp_mat[3*(j+randg)+0].x_[1];// + sftv[1];
                                        tmp_mat[3*j+0].x_[2] = tmp_mat[3*(j+randg)+0].x_[2];// + sftv[2];
                                        tmp_mat[3*j+1].x_[0] = tmp_mat[3*(j+randg)+1].x_[0];// + sftv[0];
                                        tmp_mat[3*j+1].x_[1] = tmp_mat[3*(j+randg)+1].x_[1];// + sftv[1];
                                        tmp_mat[3*j+1].x_[2] = tmp_mat[3*(j+randg)+1].x_[2];// + sftv[2];
                                        tmp_mat[3*j+2].x_[0] = tmp_mat[3*(j+randg)+2].x_[0];// + sftv[0];
                                        tmp_mat[3*j+2].x_[1] = tmp_mat[3*(j+randg)+2].x_[1];// + sftv[1];
                                        tmp_mat[3*j+2].x_[2] = tmp_mat[3*(j+randg)+2].x_[2];// + sftv[2];
                                    }
                                    int numseq=pnum;
                                //    if((rand_b-rand_a-3)<2) continue;
                                //    if(rand_a==0) continue;
                                    vector<point3f> decstr(numseq);

                                    for(int j=0;j<pnum;j++)
                                    {
                                        decstr[j].ss2=numtoss(bb[j].sst);
                                        decstr[j].ssm=numtoss(bb[j].sst);
                                        decstr[j].x=tmp_mat[3*j+1].x_[0];
                                        decstr[j].y=tmp_mat[3*j+1].x_[1];
                                        decstr[j].z=tmp_mat[3*j+1].x_[2];
                                        decstr[j].ptn.x=tmp_mat[3*j+0].x_[0];
                                        decstr[j].ptn.y=tmp_mat[3*j+0].x_[1];
                                        decstr[j].ptn.z=tmp_mat[3*j+0].x_[2];
                                        decstr[j].ptc.x=tmp_mat[3*j+2].x_[0];
                                        decstr[j].ptc.y=tmp_mat[3*j+2].x_[1];
                                        decstr[j].ptc.z=tmp_mat[3*j+2].x_[2];
                        //                tmp_tt=tmp_tt+1;
                                    }
                                    if((rand_a2)<1||(rand_b-2)<1) continue;
                                    if((rand_b2)>(pnum-1)||(rand_a+2)>(pnum-1)) continue;
                                //    mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+2);
                                //    mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b2);
                                    mcfragsweepLMP2(decstr,numseq,rand_a2,rand_b2);  
                                    tmp_tt=0;
                                    int tmp_ttx=0;
                                    for(int j=0;j<pnum;j++)
                                    {
                                        tmp_mat[3*j+1].x_[0]=decstr[j].x;
                                        tmp_mat[3*j+1].x_[1]=decstr[j].y;
                                        tmp_mat[3*j+1].x_[2]=decstr[j].z;
                                        tmp_mat[3*j+0].x_[0]=decstr[j].ptn.x;
                                        tmp_mat[3*j+0].x_[1]=decstr[j].ptn.y;
                                        tmp_mat[3*j+0].x_[2]=decstr[j].ptn.z;
                                        tmp_mat[3*j+2].x_[0]=decstr[j].ptc.x;
                                        tmp_mat[3*j+2].x_[1]=decstr[j].ptc.y;
                                        tmp_mat[3*j+2].x_[2]=decstr[j].ptc.z;
                                        tmp_tt=tmp_tt+1;
                                    }     
                                    vector<point3f>().swap(decstr);                
                                } else // right shift
                                {
                                    if((rand_a-1)<1) continue;
        //                            vector<float> sftv(3,0.0);
        //                            sftv[0]= tmp_mat[3*rand_b+1].x_[0] - tmp_mat[3*(rand_b-randg)+1].x_[0];
        //                            sftv[1]= tmp_mat[3*rand_b+1].x_[1] - tmp_mat[3*(rand_b-randg)+1].x_[1];
        //                            sftv[2]= tmp_mat[3*rand_b+1].x_[2] - tmp_mat[3*(rand_b-randg)+1].x_[2];                            
                                    for(int j=rand_b;j>=rand_a;j--)
                                    {
        //                                vector<float> sftv(3,0.0);
        //                                sftv[0]= tmp_mat[3*j+1].x_[0] - tmp_mat[3*(j-randg)+1].x_[0];
        //                                sftv[1]= tmp_mat[3*j+1].x_[1] - tmp_mat[3*(j-randg)+1].x_[1];
        //                                sftv[2]= tmp_mat[3*j+1].x_[2] - tmp_mat[3*(j-randg)+1].x_[2];
                                        tmp_mat[3*j+0].x_[0] = tmp_mat[3*(j-randg)+0].x_[0];// + sftv[0];
                                        tmp_mat[3*j+0].x_[1] = tmp_mat[3*(j-randg)+0].x_[1];// + sftv[1];
                                        tmp_mat[3*j+0].x_[2] = tmp_mat[3*(j-randg)+0].x_[2];// + sftv[2];
                                        tmp_mat[3*j+1].x_[0] = tmp_mat[3*(j-randg)+1].x_[0];// + sftv[0];
                                        tmp_mat[3*j+1].x_[1] = tmp_mat[3*(j-randg)+1].x_[1];// + sftv[1];
                                        tmp_mat[3*j+1].x_[2] = tmp_mat[3*(j-randg)+1].x_[2];// + sftv[2];
                                        tmp_mat[3*j+2].x_[0] = tmp_mat[3*(j-randg)+2].x_[0];// + sftv[0];
                                        tmp_mat[3*j+2].x_[1] = tmp_mat[3*(j-randg)+2].x_[1];// + sftv[1];
                                        tmp_mat[3*j+2].x_[2] = tmp_mat[3*(j-randg)+2].x_[2];// + sftv[2];
                                    }
                                    int numseq=pnum;
                                //    if((rand_b-rand_a-3)<2) continue;
                                //    if(rand_a==0) continue;
                                    vector<point3f> decstr(numseq);

                                    for(int j=0;j<pnum;j++)
                                    {
                                        decstr[j].ss2=numtoss(bb[j].sst);
                                        decstr[j].ssm=numtoss(bb[j].sst);
                                        decstr[j].x=tmp_mat[3*j+1].x_[0];
                                        decstr[j].y=tmp_mat[3*j+1].x_[1];
                                        decstr[j].z=tmp_mat[3*j+1].x_[2];
                                        decstr[j].ptn.x=tmp_mat[3*j+0].x_[0];
                                        decstr[j].ptn.y=tmp_mat[3*j+0].x_[1];
                                        decstr[j].ptn.z=tmp_mat[3*j+0].x_[2];
                                        decstr[j].ptc.x=tmp_mat[3*j+2].x_[0];
                                        decstr[j].ptc.y=tmp_mat[3*j+2].x_[1];
                                        decstr[j].ptc.z=tmp_mat[3*j+2].x_[2];
                        //                tmp_tt=tmp_tt+1;
                                    }
                                    if((rand_a2)<1||(rand_b-2)<1) continue;
                                    if((rand_b2)>(pnum-1)||(rand_a+2)>(pnum-1)) continue;
                        //            if((rand_a-b)<1) continue;
                        //            if((rand_b+2)>(pnum-1)) continue;
                        //            mcfragsweepLMP2(decstr,numseq,rand_a-2,rand_b+2);
                        //            mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b+2); 
                                    mcfragsweepLMP2(decstr,numseq,rand_a2,rand_b2); 
                                    tmp_tt=0;
                                    int tmp_ttx=0;
                                    for(int j=0;j<pnum;j++)
                                    {
                                        tmp_mat[3*j+1].x_[0]=decstr[j].x;
                                        tmp_mat[3*j+1].x_[1]=decstr[j].y;
                                        tmp_mat[3*j+1].x_[2]=decstr[j].z;
                                        tmp_mat[3*j+0].x_[0]=decstr[j].ptn.x;
                                        tmp_mat[3*j+0].x_[1]=decstr[j].ptn.y;
                                        tmp_mat[3*j+0].x_[2]=decstr[j].ptn.z;
                                        tmp_mat[3*j+2].x_[0]=decstr[j].ptc.x;
                                        tmp_mat[3*j+2].x_[1]=decstr[j].ptc.y;
                                        tmp_mat[3*j+2].x_[2]=decstr[j].ptc.z;
                                        tmp_tt=tmp_tt+1;
                                    }     
                                    vector<point3f>().swap(decstr);                                        
                                }
                            } 
                    //        if(rand_movement==7)
                    //        {
//                                theDensityMap.select_points( tmp_mat, theDensityMap.density);
//                                int nnum = theDensityMap.points_to_search_.size();                                  
                    //            for(int j=0;j<nnum;i++)
                    //        }                                                                 
                        }
        //                cout<<" 1 ";
        /*                int randxx = rand()%(rand_b-rand_a+1)+rand_a;
                        vector<float> tmpcc(3,0.0);
                        int randyy=randxx-rand_a;  // randxx and randyy is the same point
                        float deta0 = (RandomDoubleX(0.0,1.0)*2.0-1.0)*0.1+0.0;
                        float deta1 = (RandomDoubleX(0.0,1.0)*2.0-1.0)*0.1+0.0;
                        float deta2 = (RandomDoubleX(0.0,1.0)*2.0-1.0)*0.1+0.0;
                        tmpcc[0] = tmp_mat[randyy].x_[0] + deta0;
                        tmpcc[1] = tmp_mat[randyy].x_[1] + deta1;
                        tmpcc[2] = tmp_mat[randyy].x_[2] + deta2;
                        float dis_s0=0.0,dis_s1=0.0;
                        if(randxx==0) {
                            dis_s0 = Distance_point(tmpcc,pointsBx[randxx+1].x_);
                        }
                        else if(randxx==(pnum-1))
                        {
                            dis_s0 = Distance_point(tmpcc,pointsBx[randxx-1].x_);
                        }
                        else {
                            if(randyy-1>=0) dis_s0 = Distance_point(tmpcc,tmp_mat[randyy-1].x_);
                            if(randyy+1 < tmp_mat.size()) dis_s0 = Distance_point(tmpcc,tmp_mat[randyy+1].x_);
                            if(randyy==(tmp_mat.size()-1)&& (randxx+1)<=(pnumx-1)) dis_s0 = Distance_point(tmpcc,pointsBx[randxx+1].x_);
                            if(randyy==0&&(randxx-1)>=0) dis_s0 = Distance_point(tmpcc,pointsBx[randxx-1].x_);
                        } */
                //      if(dis_s0>3.70 && dis_s0<3.85 ) tmp_mat[randyy].x_ = tmpcc;   
                //        if(dis_s0>1.43001 && dis_s0<1.49999 && dis_s1<(1.33999) && dis_s1>(1.31001)) tmp_mat[randyy] = tmpcc;
        //              cout<<" 2 ";

                        dis_p=0.0;
                        tp_clashx= false;
                        new_vwd = 0.0;
                        new_clash = 0.0;
                        for(int jj=3*rand_a2;jj<=3*rand_b2+2;jj++)
                        {
                            if(tmp_mat[jj].elt_=="CA")
                            {
                    //            float VR1= GetVdwRadius(tmp_mat[jj].elt_);
                                string VR1 = tmp_mat[jj].elt_;
                            //    for(int t=0;t<pointsBx.size();t++)
                                for(int t=3*dm_s+0;t<3*dm_e;t++)
                                {
                                    if(tmp_mat[t].elt_=="CA")
                                    {
                        //                if((jj+rand_a)==t) continue;
                                        if(abs(jj-t)<3) continue;
                        //                if(t>3*rand_a && t<(3*rand_a+tmp_mat.size())) continue;
                        //                float VR2= GetVdwRadius(pointsBx[t].elt_);
                                        string VR2 = tmp_mat[t].elt_ ;
                                        dis_p = Distance_point(tmp_mat[jj].x_,tmp_mat[t].x_);
                                        //if(dis_p < (VR1+VR2))
                                //        if(dis_p < 3.70)
                                //        {
                         //                   new_clash = new_clash + 1.0/ sqrt(dis_p+0.0001) ;
                         //                   old_clash = old_clash + exp( 3.75-dis_p) ;
                                            new_clash = new_clash + GetVdwEgC(VR1,VR2,dis_p);
                                //        }                        
                        //                new_vwd = new_vwd + GetVdwEg(VR1,VR2,dis_p);
                        //                new_vwd = new_vwd + GetVdwEgC(VR1,VR2,dis_p);
                        //                if((dis_p)<2.0) tp_clashx = true;
                        //                if(tp_clashx)
                        //                {
                        //                    cout<<"jj,t: "<<jj+rand_a<<" "<<t<<" dis: "<<VR1+VR2<<" "<<dis_p<<endl;
                        //                    tp_clashx=false;
                        //                }
                    //                    cout<< "VR1+VR2,dis_p: "<<VR1+VR2<<" "<<dis_p<<endl;
                                    }
                                }
                //                if(tp_clashx) break;
                            }
                        }
            //            cout<<"tp_clashx: " <<tp_clashx<<endl;
            //            if(tp_clashx) continue;  
            //            cout<<" 3 "<<endl;                 

                    for ( int t=0; t<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++t ) {
                        rhoC2[t]=0.0;
                        inv_rho_mask2[t]=1.0;
                    }
                    del_ij=vector<float>(3,0.0);
                    atm_j=vector<float>(3,0.0);
                    tm_i1 =0;
                    atm_idy=vector<float>(3,0.0);
                    bond0max=-1000.0;bond1max=-1000.0;bond2max=-1000.0;
                    bond0min=1000.0;bond1min=1000.0;bond2min=1000.0;

                    for(int j=3*rand_a2;j<=3*rand_b2+2;j++)
                    {
                        if(tmp_mat[j].elt_ == "CA")
                        {
                            vector<float> cartX1;
                            vector<float> fracX1; 
                        //    elt_i = tmp_mat[3*tm_i1+1].elt_;
                            elt_i = tmp_mat[j].elt_;
                            elt_i = elt_i[0];
                            OneGaussianScattering sig_j = get_A( elt_i );
                            k = sig_j.k( theDensityMap.effectiveB );
                            C = sig_j.C( k );
                            if ( C < 1e-6 ) continue;  

                    //        cartX1 = tmp_mat[3*tm_i1+1].x_;
                            cartX1 = tmp_mat[j].x_;
                    //        tm_i1 = tm_i1 + 1;
                            MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);
                            atm_idy[0] = pos_mod (double(fracX1[0]*grid[0] - origin[0] + 1) , (double)grid[0]);
                            atm_idy[1] = pos_mod (double(fracX1[1]*grid[1] - origin[1] + 1) , (double)grid[1]);
                            atm_idy[2] = pos_mod (double(fracX1[2]*grid[2] - origin[2] + 1) , (double)grid[2]);                 
                            for(int z=1;z<=theDensityMap.density.u3();z++)
                            {
                                atm_j[2] = z;
                                del_ij[2] =(atm_idy[2]-atm_j[2])/grid[2];
                                if ( del_ij[2] > 0.5 ) del_ij[2]-=1.0;
                                if ( del_ij[2] < -0.5 ) del_ij[2]+=1.0;                    
                                del_ij[0] = del_ij[1] = 0.0;
                                vector<float> frac_tmpz;
                                MatrixTimesTransVector(theDensityMap.f2c,del_ij,frac_tmpz);
                                if(square_len(frac_tmpz)> (ppadding+atom_m)*(ppadding+atom_m) ) continue;
                                if(z<bond2min) bond2min = z; 
                                if(z>bond2max) bond2max = z;              
                                for(int y=1;y<=theDensityMap.density.u2();y++)
                                {
                                    atm_j[1] = y;
                                    del_ij[1] = (atm_idy[1] - atm_j[1])/grid[1];
                                    // wrap-around??
                                    if ( del_ij[1] > 0.5 ) del_ij[1]-=1.0;
                                    if ( del_ij[1] < -0.5 ) del_ij[1]+=1.0;                        
                                    del_ij[0] = 0.0;
                                    vector<float> frac_tmpy;
                                    MatrixTimesTransVector(theDensityMap.f2c,del_ij,frac_tmpy);
                                    if(square_len(frac_tmpy)> (ppadding+atom_m)*(ppadding+atom_m) ) continue; 
                                    if(y>bond1max) bond1max = y; 
                                    if(y<bond1min) bond1min = y;                    
                                    for(int x=1;x<=theDensityMap.density.u1();x++)
                                    {
                                        atm_j[0] = x;
                                        del_ij[0] = (atm_idy[0] - atm_j[0])/grid[0];
                                        // wrap-around??
                                        if ( del_ij[0] > 0.5 ) del_ij[0]-=1.0;
                                        if ( del_ij[0] < -0.5 ) del_ij[0]+=1.0;                            
                                        vector<float> cart_del_ij2;
                                        MatrixTimesTransVector(theDensityMap.f2c,del_ij,cart_del_ij2);
                                        float d2 = square_len(cart_del_ij2);
                                        if(d2 > (ppadding+atom_m)*(ppadding+atom_m) ) continue;
                                    
                                        float atm = C*exp(-k*d2);
                                        float sigmoid_msk = exp( d2 - (atom_m)*(atom_m)  );
                                    //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK));
                                    //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK));
                                        float inv_msk = 1.0/(1.0+sigmoid_msk);
                        //                rhoC2(x,y,z) += atm;
                        //                rhoC2(x,y,z) += atm;
                                        rhoC2(x,y,z) += atm;
                                        inv_rho_mask2(x,y,z) *= (1.0 - inv_msk);


                                        if(x>bond0max) bond0max = x;
                                        if(x<bond0min) bond0min = x;
                                        
                                        
            //                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
            //                                mask2(x,y,z) = 1.0; // problem?
            //                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
            //                                    float atm = C*exp(-k*d2);
            //                                    rhoC(x,y,z) += atm;
            //                                }
            //                            }
                                    }
                                }
                            }
                            tm_i1 = tm_i1 + 1;
                        }
                    }
            //        rhoC2 = rhoC01;
                    max_rhc = 0.0;
        //            for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
        //                if(rhoC2[x]>max_rhc) max_rhc = rhoC2[x] ;
        //            }     
            /*        for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
                //        cout<<theDensityMap.density[x]<<" ";
                        if(rhoC2[x]<(3.0/5.0)*max_rhc) rhoC2[x]=0; 
                    } */
        //            rhoC2 = rhoC01;                           
                    sumC_i2=0.0, sumO_i2=0.0, sumCO_i2=0.0, vol_i2=0.0, CC_i2=0.0;
                    sumO2_i2=0.0, sumC2_i2=0.0, varC_i2=0.0, varO_i2=0.0;
                    clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
                /*        for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
                //            tmp_den=0.0;
                //            tmp_rhc =0.0;
                //            if(theDensityMap.density[x] > (2.0/3.0)*max_den) tmp_den= theDensityMap.density[x];
                //            if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
                //            if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                            // fetch this point
                            clc_x2 = rhoC2[x];
                //            clc_x2 = tmp_rhc;
                            obs_x2 = theDensityMap.density[x];
                //            obs_x2 = tmp_den;
                            eps_x2 = 1-inv_rho_mask2[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                //            eps_x2 = 1.0;

                            // SMOOTHED
                            sumCO_i2 += eps_x2*clc_x2*obs_x2;
                            sumO_i2  += eps_x2*obs_x2;
                            sumO2_i2 += eps_x2*obs_x2*obs_x2;
                            sumC_i2  += eps_x2*clc_x2;
                            sumC2_i2 += eps_x2*clc_x2*clc_x2;
                            vol_i2   += eps_x2;
                        }
                        varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
                        varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
                        if ( varC_i2 == 0 || varO_i2 == 0 ) {
                            CC_i2 = 0;
                        } else {
                            CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
                        }   */
                    for(int x=bond0min;x<=bond0max;x++)
                    {
                       for(int y=bond1min;y<=bond1max;y++)
                        {
                            for(int z=bond2min;z<=bond2max;z++)
                            {
                   //            tmp_den=0.0;
                    //            tmp_rhc =0.0;
                    //            if(theDensityMap.density[x] > (2.0/3.0)*max_den) tmp_den= theDensityMap.density[x];
                    //            if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
                    //            if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                                // fetch this point
                                clc_x2 = rhoC2(x,y,z);
                    //            clc_x2 = tmp_rhc;
                                obs_x2 = theDensityMap.density(x,y,z);
                    //            obs_x2 = tmp_den;
                                eps_x2 = 1-inv_rho_mask2(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                    //            eps_x2 = 1.0;

                                // SMOOTHED
                                sumCO_i2 += eps_x2*clc_x2*obs_x2;
                                sumO_i2  += eps_x2*obs_x2;
                                sumO2_i2 += eps_x2*obs_x2*obs_x2;
                                sumC_i2  += eps_x2*clc_x2;
                                sumC2_i2 += eps_x2*clc_x2*clc_x2;
                                vol_i2   += eps_x2;
                            }
                        }
                    }          
                    varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
                    varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
                    if ( varC_i2 == 0 || varO_i2 == 0 ) {
                        CC_i2 = 0;
                    } else {
                        CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
                    }       
                    new_CC =1.0000-CC_i2;       
            /*        sumis=0.0;
                    sumu=0.0;
                    for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
                        tmp_den = 0.0; 
                        tmp_rhc=0.0;
                        if(theDensityMap.density[x]>(3.0/5.0)*max_den) tmp_den=1.0;
                        if(rhoC2[x]>(3.0/5.0)*max_rhc) tmp_rhc=1.0;
        //                if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                        // fetch this point
                        sumis=sumis+tmp_den*tmp_rhc;
                        sumu=sumu+tmp_den+tmp_rhc-tmp_den*tmp_rhc;
                    }
                    if ( sumu == 0 ) {
                        CC_i2 = 0;
                    } else {
                        CC_i2 = sumis/sumu;
                    }   */
            /*        sumis=0.0;
                    sumu=0.0;
                    for(int x=bond0min;x<=bond0max;x++)
                    {
                       for(int y=bond1min;y<=bond1max;y++)
                        {
                            for(int z=bond2min;z<=bond2max;z++)
                            {
                                tmp_den = 0.0; 
                                tmp_rhc=0.0;
                                if(theDensityMap.density(x,y,z)>(3.0/5.0)*max_den) tmp_den=1.0;
                                if(rhoC2(x,y,z)>(3.0/5.0)*max_rhc) tmp_rhc=1.0;
                //                if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                                // fetch this point
                                sumis=sumis+tmp_den*tmp_rhc;
                                sumu=sumu+tmp_den+tmp_rhc-tmp_den*tmp_rhc;
                            }
                        }
                    }          
                    if ( sumu == 0 ) {
                        CC_i2 = 0;
                    } else {
                        CC_i2 = sumis/sumu;
                    }             */
        /*            sumis=0.0;
                    sumu=0.0;   
                    bond0max=0.0;bond1max=0.0;bond2max=0.0;
                    bond0min=1000.0;bond1min=1000.0;bond2min=1000.0;
                    tm_i1=0;
                    atm_idy=vector<vector<float> > ((rand_b-rand_a+1),vector<float>(3,0.0));
                    for(int j=0;j<tmp_mat.size();j++)
                    {
         //               tmp_mat.push_back(pointsBx[j]);
                        if(tmp_mat[j].elt_ == "CA")
                        {
                            vector<float> cartX1;
                            vector<float> fracX1; 

                     //       cartX1 = tmp_mat[3*tm_i1+1].x_;
                            cartX1 = tmp_mat[j].x_;
                            
                            MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);
                            atm_idy[tm_i1][0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                            atm_idy[tm_i1][1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                            atm_idy[tm_i1][2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);
                            if(atm_idy[tm_i1][0]>bond0max) bond0max = atm_idy[tm_i1][0];
                            if(atm_idy[tm_i1][1]>bond1max) bond1max = atm_idy[tm_i1][1];
                            if(atm_idy[tm_i1][2]>bond2max) bond2max = atm_idy[tm_i1][2];
                            if(atm_idy[tm_i1][0]<bond0min) bond0min = atm_idy[tm_i1][0];
                            if(atm_idy[tm_i1][1]<bond1min) bond1min = atm_idy[tm_i1][1];
                            if(atm_idy[tm_i1][2]<bond2min) bond2min = atm_idy[tm_i1][2];
                            tm_i1 = tm_i1 + 1;
                        }
                    }  
                    for(int x=bond0min;x<=bond0max;x++)
                    {
                       for(int y=bond1min;y<=bond1max;y++)
                        {
                            for(int z=bond2min;z<=bond2max;z++)
                            {
                                tmp_den = 0.0; 
                                tmp_rhc=0.0;
                                if(theDensityMap.density(x,y,z)>(3.0/5.0)*max_den) tmp_den=1.0;
                                if(rhoC2(x,y,z)>(3.0/5.0)*max_rhc) tmp_rhc=1.0;
                //                if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                                // fetch this point
                                sumis=sumis+tmp_den*tmp_rhc;
                                sumu=sumu+tmp_den+tmp_rhc-tmp_den*tmp_rhc;
                            }
                        }
                    }          
                    if ( sumu == 0 ) {
                        CC_i2 = 0;
                    } else {
                        CC_i2 = sumis/sumu;
                    }      */                         
        //            vector<vector<float> > bondx(3,vector<float>(2,0.0));
        /*            bond0max=0.0;bond1max=0.0;bond2max=0.0;
                    bond0min=1000.0;bond1min=1000.0;bond2min=1000.0;
                    tm_i1=0;
                    atm_idy=vector<vector<float> > ((rand_b-rand_a+1),vector<float>(3,0.0));
                    for(int j=0;j<tmp_mat.size();j++)
                    {
         //               tmp_mat.push_back(pointsBx[j]);
                        if(tmp_mat[j].elt_ == "CA")
                        {
                            vector<float> cartX1;
                            vector<float> fracX1; 

                     //       cartX1 = tmp_mat[3*tm_i1+1].x_;
                            cartX1 = tmp_mat[j].x_;
                            
                            MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);
                            atm_idy[tm_i1][0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                            atm_idy[tm_i1][1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                            atm_idy[tm_i1][2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);
                            if(atm_idy[tm_i1][0]>bond0max) bond0max = atm_idy[tm_i1][0];
                            if(atm_idy[tm_i1][1]>bond1max) bond1max = atm_idy[tm_i1][1];
                            if(atm_idy[tm_i1][2]>bond2max) bond2max = atm_idy[tm_i1][2];
                            if(atm_idy[tm_i1][0]<bond0min) bond0min = atm_idy[tm_i1][0];
                            if(atm_idy[tm_i1][1]<bond1min) bond1min = atm_idy[tm_i1][1];
                            if(atm_idy[tm_i1][2]<bond2min) bond2min = atm_idy[tm_i1][2];
                            tm_i1 = tm_i1 + 1;
                        }
                    } 
                    for(int x=bond0min;x<=bond0max;x++)
                    {
                       for(int y=bond1min;y<=bond1max;y++)
                        {
                            for(int z=bond2min;z<=bond2max;z++)
                            {
                                clc_x2 = rhoC2(x,y,z);
                //                clc_x2 = tmp_rhc;
                                obs_x2 = theDensityMap.density(x,y,z);
                //                obs_x2 = tmp_den;
                //                eps_x2 = 1-inv_rho_mask2[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                                eps_x2 = 1.0;

                                // SMOOTHED
                                sumCO_i2 += eps_x2*clc_x2*obs_x2;
                                sumO_i2  += eps_x2*obs_x2;
                                sumO2_i2 += eps_x2*obs_x2*obs_x2;
                                sumC_i2  += eps_x2*clc_x2;
                                sumC2_i2 += eps_x2*clc_x2*clc_x2;
                                vol_i2   += eps_x2;
                            }
                        }
                    }          
                    varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
                    varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
                    if ( varC_i2 == 0 || varO_i2 == 0 ) {
                        CC_i2 = 0;
                    } else {
                        CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
                    }                */
         /*               for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
                            clc_x2 = rhoC2[x];
                            obs_x2 = theDensityMap.density[x];
                            eps_x2 = 1-inv_rho_mask2[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal

                            // SMOOTHED
                            sumCO_i2 += eps_x2*clc_x2*obs_x2;
                            sumO_i2  += eps_x2*obs_x2;
                            sumO2_i2 += eps_x2*obs_x2*obs_x2;
                            sumC_i2  += eps_x2*clc_x2;
                            sumC2_i2 += eps_x2*clc_x2*clc_x2;
                            vol_i2   += eps_x2;
                        }
                        sumO_i2 = sumO_i2/vol_i2;
                        sumC_i2 = sumC_i2/vol_i2;
                        sumCO_i2 = sumCO_i2/vol_i2;
                        sumC2_i2 = sumC2_i2/vol_i2;
                        sumO2_i2 = sumO2_i2/vol_i2;
                        varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2);
                        varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2) ;
                        if ( varC_i2 == 0 || varO_i2 == 0 ) {
                            CC_i2 = 0;
                        } else {
                            CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2) / sqrt( varC_i2 * varO_i2 );
                        }  */
            /*        for(int x=bond0min;x<=bond0max;x++)
                    {
                       for(int y=bond1min;y<=bond1max;y++)
                        {
                            for(int z=bond2min;z<=bond2max;z++)
                            {
                                clc_x2 = rhoC2(x,y,z);
                //                clc_x2 = tmp_rhc;
                                obs_x2 = theDensityMap.density(x,y,z);
                //                obs_x2 = tmp_den;
                //                eps_x2 = 1-inv_rho_mask2[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                                eps_x2 = 1.0;

                                // SMOOTHED
                                vol_i2 += eps_x2;
                                sumC_i2 += (eps_x2*clc_x2*clc_x2);
                                sumO_i2 += (eps_x2*obs_x2*obs_x2);
                                sumCO_i2 += (eps_x2*clc_x2*obs_x2);
                            }
                        }
                    }          
                    if ( sumC_i2 == 0 || sumO_i2 == 0 ) {
                        CC_i2 = 0;
                    } else {
                        CC_i2 = sumCO_i2/(sqrt(sumC_i2+sumO_i2)); 
                    }             */
        /*                for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x )
                        {
                            eps_x2 = 1-inv_rho_mask2[x];
                            vol_i2 += eps_x2;
                            sumC_i2 += (eps_x2*rhoC2[x]*rhoC2[x]);
                            sumO_i2 += (eps_x2*theDensityMap.density[x]*theDensityMap.density[x]);
                            sumCO_i2 += (eps_x2*rhoC2[x]*theDensityMap.density[x]);
                        }
                    //    sumCO_i2 = sumCO_i2/vol_i2;
                    //    sumO_i2 = sumO_i2/vol_i2;
                    //    sumC_i2 = sumC_i2/vol_i2;
                        CC_i2 = sumCO_i2/(sqrt(sumC_i2*sumO_i2));  */

        //                cout<< "DDDDDXX"<<endl;
                        // ss inforamtion     
                        float new_Ehbond=0.0;
                //        point3f *new_decstr;

                    /*    vector<boneinfo> bby(numseq);
                        tmp_tt=0;
                        for(int j=rand_a;j<=rand_b;j++)
                        {
                            bby[tmp_tt]=bb[j];
                            tmp_tt = tmp_tt +1;
                        }       */  
                     //   if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E')
                        {                       
                    /*        tmp_tt=0;
                            for(int j=0;j<pnum;j++)
                            {
                                bby[j].indn=3*j+0;
                                bby[j].indca=3*j+1;
                                bby[j].indc=3*j+2; 
                                string resng = pointsrenm[j];
                                char resny = resng[0];
                                bby[j].resid = aminoid(resny);                                    
                                tmp_tt=tmp_tt+1;
                            }     
                        //    calcsse2(bby, pnum, tmp_mat); 
                            calcssennhoc(bby,pnum,tmp_mat);   */
                            vector<boneinfo> bbg(pnum);
                            for(int j=0;j<pnum;j++)
                            {
                                bbg[j].indn = 3*j+0;
                                bbg[j].indca = 3*j+1;
                                bbg[j].indc = 3*j+2; 
                                string resng = pointsrenm[j];
                                char resny = resng[0];
                                bbg[j].resid = aminoid(resny);                                    
                            }     
                        //    calcsse2(bbg, pnum, tmp_mat); 
                            calcssennhoc(bbg,pnum,tmp_mat);                

        //                if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E')
        //                {
                            int numseq=pnum;
                            vector<point3f> new_decstr(numseq);             
                            tmp_tt=0;
                            for(int j=0;j<pnum;j++)
                            {                     
                            //    new_decstr[tmp_tt].ss2=numtoss(bby[j].sst);
                                new_decstr[j].ss2=numtoss(bbg[j].sst);
                                new_decstr[j].ssm=numtoss(bbg[j].sst);
                                new_decstr[j].x=tmp_mat[3*tmp_tt+1].x_[0];
                                new_decstr[j].y=tmp_mat[3*tmp_tt+1].x_[1];
                                new_decstr[j].z=tmp_mat[3*tmp_tt+1].x_[2];
                                new_decstr[j].ptn.x=tmp_mat[3*tmp_tt+0].x_[0];
                                new_decstr[j].ptn.y=tmp_mat[3*tmp_tt+0].x_[1];
                                new_decstr[j].ptn.z=tmp_mat[3*tmp_tt+0].x_[2];
                                new_decstr[j].ptc.x=tmp_mat[3*tmp_tt+2].x_[0];
                                new_decstr[j].ptc.y=tmp_mat[3*tmp_tt+2].x_[1];
                                new_decstr[j].ptc.z=tmp_mat[3*tmp_tt+2].x_[2];
                                tmp_tt=tmp_tt+1;
                            }                  
                        //    new_Ehbond=energyhbondcanc(new_decstr,numseq);
                            new_Ehbond = energyhbondnhoc2(new_decstr,numseq);                           
                            vector<point3f>().swap(new_decstr);
                        }        
                    //    cout<<"new_Ehbond: "<<new_Ehbond<<endl; 
                        tm_i1=0;                    
            /*            float new_ss_E =0.0;
                        vector<string > new_ss_inf(CA_num);
                        vector<vector<float> > new_CA_cor;
                        tm_i1=0;                
                        for(int j=0;j<tmp_mat.size();j++)
                        {
                            if(tmp_mat[j].elt_ == "CA")
                            {
                                new_CA_cor.push_back(tmp_mat[j].x_);
                                tm_i1 = tm_i1+1 ; 
                            }                 
                        }
                        make_secx(new_CA_cor,CA_num,new_ss_inf);
                        ss_ab_num=0; // helix and b strand number
                        ss_c_num =0; // coil number                  
                        for(int j=0;j<new_ss_inf.size();j++)
                        {
                            if(new_ss_inf[j] == 'H')
                            {
                                new_ss_E = new_ss_E + 0.1;
                                ss_ab_num = ss_ab_num +1;
                            }
                            else if(new_ss_inf[j] == 'E')
                            {
                                new_ss_E = new_ss_E + 0.1;
                                ss_ab_num = ss_ab_num +1;
                            }
                            else{
                                new_ss_E = new_ss_E + 0.2;
                                ss_c_num = ss_c_num +1;
                            }                    
                        }
                        new_ss_E = new_ss_E/(CA_num);    */             
                        
                /*        vector<point3f> decstrbondangx(pnum);
                        for(int j=0;j<pnum;j++)
                        {
                            decstrbondangx[j].ss2=numtoss(bb[j].sst);
                            decstrbondangx[j].x=tmp_mat[3*j+1].x_[0];
                            decstrbondangx[j].y=tmp_mat[3*j+1].x_[1];
                            decstrbondangx[j].z=tmp_mat[3*j+1].x_[2];
                            decstrbondangx[j].ptn.x=tmp_mat[3*j+0].x_[0];
                            decstrbondangx[j].ptn.y=tmp_mat[3*j+0].x_[1];
                            decstrbondangx[j].ptn.z=tmp_mat[3*j+0].x_[2];
                            decstrbondangx[j].ptc.x=tmp_mat[3*j+2].x_[0];
                            decstrbondangx[j].ptc.y=tmp_mat[3*j+2].x_[1];
                            decstrbondangx[j].ptc.z=tmp_mat[3*j+2].x_[2];
                            string resng = pointsrenm[j];
                            char resny = resng[0];
                            decstrbondangx[j].iaa = aminoid(resny);
                //                decstr[i].iaa=pp.aminoid(decstr[i].aaa);                
                //                tmp_tt=tmp_tt+1;
                        }     
                        for(int j=0;j<pnum;j++)
                        {
                            if(j==0)
                            {
                                decstrbondangx[j].ang[0]=0.0;
                                decstrbondangx[j].ang[1]=0.0;
                                decstrbondangx[j].ang[2]=0.0;
                            }               
                            point3d p12 = setv(decstrbondangx[j-1].x-decstrbondangx[j-1].ptc.x,decstrbondangx[j-1].y-decstrbondangx[j-1].ptc.y,decstrbondangx[j-1].z-decstrbondangx[j-1].ptc.z);
                            point3d p23 = setv(decstrbondangx[j].ptn.x-decstrbondangx[j-1].ptc.x,decstrbondangx[j].ptn.y-decstrbondangx[j-1].ptc.y,decstrbondangx[j].ptn.z-decstrbondangx[j-1].ptc.z);
                            decstrbondangx[j].ang[0]=angv(p12,p23)*degrad;
                            p12 = setv(decstrbondangx[j-1].ptc.x-decstrbondangx[j].ptn.x,decstrbondangx[j-1].ptc.y-decstrbondangx[j].ptn.y,decstrbondangx[j-1].ptc.z-decstrbondangx[j].ptn.z);
                            p23 = setv(decstrbondangx[j].x-decstrbondangx[j].ptn.x,decstrbondangx[j].y-decstrbondangx[j].ptn.y,decstrbondangx[j].z-decstrbondangx[j].ptn.z);
                            decstrbondangx[j].ang[1]=angv(p12,p23)*degrad;
                            p12 = setv(decstrbondangx[j].ptn.x-decstrbondangx[j].x,decstrbondangx[j].ptn.y-decstrbondangx[j].y,decstrbondangx[j].ptn.z-decstrbondangx[j].z);
                            p23 = setv(decstrbondangx[j].ptc.x-decstrbondangx[j].x,decstrbondangx[j].ptc.y-decstrbondangx[j].y,decstrbondangx[j].ptc.z-decstrbondangx[j].z);
                            decstrbondangx[j].ang[2]= angv(p12,p23)*degrad;
                //                tmp_tt=tmp_tt+1;
                        }
                        float new_bondangenergy =0.0;
                        Eenergybondangle(decstrbondangx,pnum,new_bondangenergy);
                        cout<< "new_bondangenergy: "<< new_bondangenergy/pnum<<endl;   */

                //        cout<< "new_CC: "<<new_CC<<" ";
                    //    new_dE = 30*new_CC + new_vwd/((pnumx-1)*(rand_b-rand_a+1));
                    //    new_dE = 0.8*new_CC + 0.5*new_clash/abs(rand_b-rand_a+1) + 0.5*new_Ehbond/abs(rand_b-rand_a+1);
                    //    new_dE = 500.0*new_CC + 100.0*new_clash/abs(rand_b2-rand_a2+1) + 100.0*(new_Ehbond/abs(pnum))/(1.0-(new_Ehbond/abs(pnum)));// + new_bondangenergy/pnum;
                        float Eres=0.0;
                   /*     for(int j=0;j<pnum;j++)
                        {
                            Eres = Eres + Distance_point(tmp_mat[3*j+1].x_,pointsBy[3*j+1].x_);
                        }. */
                        for(int j=3*rand_a2;j<=3*rand_b2+2;j++)
                        {
                            for(int js=j-3;js<=j+3;js++)
                            {
                                if(js==j) continue;
                                if(js<0) continue;
                                if(js>(pnum-1)) continue;
                                float dx0 = Distance_point(tmp_mat[j].x_,tmp_mat[js].x_);
                                float dx1 = Distance_point(pointsBx0[j].x_,pointsBx0[js].x_);
                                Eres = Eres + (dx0-dx1)*(dx0-dx1);
                            }
                        }                        
                        Eres = Eres/float(rand_b2-rand_a2+1); 

                        vector<point3f>().swap(decstrbondlen);
                        decstrbondlen = vector<point3f>(pnum);
                        for(int j=0;j<pnum;j++)
                        {
                            decstrbondlen[j].ss2=numtoss(bb[j].sst);
                            decstrbondlen[j].ssm=numtoss(bb[j].sst);
                            decstrbondlen[j].stype=numtoss(bb[j].sst);
                            decstrbondlen[j].x=tmp_mat[3*j+1].x_[0];
                            decstrbondlen[j].y=tmp_mat[3*j+1].x_[1];
                            decstrbondlen[j].z=tmp_mat[3*j+1].x_[2];
                            decstrbondlen[j].ptn.x=tmp_mat[3*j+0].x_[0];
                            decstrbondlen[j].ptn.y=tmp_mat[3*j+0].x_[1];
                            decstrbondlen[j].ptn.z=tmp_mat[3*j+0].x_[2];
                            decstrbondlen[j].ptc.x=tmp_mat[3*j+2].x_[0];
                            decstrbondlen[j].ptc.y=tmp_mat[3*j+2].x_[1];
                            decstrbondlen[j].ptc.z=tmp_mat[3*j+2].x_[2];
                            string resng = pointsrenm[j];
                            char resny = resng[0];
                            decstrbondlen[j].iaa = aminoid(resny);
            //                decstr[i].iaa=pp.aminoid(decstr[i].aaa);                
            //                tmp_tt=tmp_tt+1;
                        } 
                        str2tor(decstrbondlen,pnum,3);
                    //    str2torp(decstrbondlen,pnum,0,pnum-1); 
                /*        float new_fang=0.0;
                        float new_Ebondang=0.0;
                        new_Ebondang =(float) Eenergybondangle(decstrbondlen,pnum,new_fang);
                        new_Ebondang = 0.30*new_Ebondang + new_fang;   */                    
                        float new_fene=0.0;
                        float new_Ebondlen=0.0;
                        new_Ebondlen = (float) energybondlength(decstrbondlen,pnum,new_fene);
                        new_Ebondlen = 0.5*new_Ebondlen + 20.0*new_fene; 
                        // torsion angle
                        float new_tor_fene=0.0;
                        float new_tor_E=0.0;
                        new_tor_E = energyrama(decstrbondlen,pnum,new_tor_fene,rlogduke,ramaduke); 
                        new_tor_E = 4.00*new_tor_E + new_tor_fene; 

                        float E_distx = 0.0;
                        for(int iu=0;iu<N_dist;iu++)
                        {
                            int ik = dist_pre[iu].ires;
                            int jk = dist_pre[iu].jres;
                            float tmp_x = dist_pre[iu].dv;
                            float tmp_y = Distance_point(tmp_mat[3*ik+1].x_,tmp_mat[3*jk+1].x_);
                    //        float tmp_z=sqrt((tmp_x-tmp_y)*(tmp_x-tmp_y))-dist_pre[iu].ave;
                            float tmp_z=tmp_y-dist_pre[iu].ave;
                            E_distx = E_distx + (log(tmp_z*tmp_z+1))/sqrt(dist_pre[iu].stdx);
                        }    
                        E_distx=E_distx/float(pnum*(rand_b2-rand_a2+1));

                        float E_conx = 0.0;
                        for(int iu=0;iu<N_distg;iu++)
                        {
                            int ik = dist_map[iu].ires;
                            int jk = dist_map[iu].jres;
                            float tmp_x = dist_map[iu].dv;
                            float tmp_y = Distance_point(tmp_mat[3*ik+1].x_,tmp_mat[3*jk+1].x_);
                            E_conx = E_conx + sqrt((tmp_x-tmp_y)*(tmp_x-tmp_y));
                        }
                        E_conx=E_conx/float((rand_b2-rand_a2+1));

                        float segdisteng=0.0;  
                        if(TF_segmx)
                        {
                            int segint=segm.size();
                            int ttu=1;
                            for(int il=0;il<segint;il++)
                            {
                        //        int intu=segnum[il];
                                int il0=segm[il].init;
                                int iln=segm[il].term;
                                for(int iul=il0;iul<iln;iul++)
                                {
                                    segdisteng=segdisteng+Distance_point(tmp_mat[3*iul+1].x_,pointsBy[3*iul+1].x_);
                                    ttu=ttu+1;
                                }
                            }
                            segdisteng=segdisteng/float(ttu);
                            segdisteng=-segdisteng;  // care wu qiong yuan
                    //        segdisteng=1.0/(segdisteng+1.0);
                        }                   
                /*        if(TF_segmx)
                        {
                            int segint=segnum.size();
                            int ttu=1;
                            for(int il=0;il<segint;il++)
                            {
                                int intu=segnum[il];
                                int il0=segm[intu].init;
                                int iln=segm[intu].term;
                                for(int iul=il0;iul<iln;iul++)
                                {
                                    segdisteng=segdisteng+Distance_point(tmp_mat[3*iul+1].x_,pointsBy[3*iul+1].x_);
                                    ttu=ttu+1;
                                }
                            }
                            segdisteng=segdisteng/float(ttu);
                            segdisteng=-segdisteng;  // care wu qiong yuan
                    //        segdisteng=1.0/(segdisteng+1.0);
                        }。*/
                        new_dE = w1*new_CC + 1.0*new_clash + 1.0*(new_Ehbond) + new_Ebondlen + 1.0*E_conx + 1.0*new_tor_E + w6*E_distx+segdisteng+Eres;// + 0.2*Eres+ new_Ebondang/500.0;

                        vector<point3f>().swap(decstrbondlen);
            //            cout<<"old_fene,new_fene: "<<old_Ebondlen<<" "<<new_Ebondlen<<endl;
                     //   new_dE = new_CC + new_clash + (new_Ehbond/abs(pnum))/(1.0-(new_Ehbond/abs(pnum)));
                    //    cout<<"new coor:  "<<tmp_mat[5].x_[0]<<" "<<tmp_mat[5].x_[1]<<" "<<tmp_mat[5].x_[2]<<endl;
                   //     cout<<"tttz: "<<200.0*old_CC<<" "<<1.0*old_clash<<" "<<10.0*(old_Ehbond)<<endl;
                    //    cout<<"ttty: "<<200.0*new_CC<<" "<<1.0*new_clash<<" "<<10.0*(new_Ehbond)<<endl;
                    //    cout<<"old: "<< old_CC<<" " <<old_clash<<" "<<old_Ehbond<<endl;
        //                cout<<"new: "<< new_CC<<" " <<new_clash<<" "<<new_Ehbond<<endl;
            //            cout<<"dE_CC,clash,Hbond,bondlen: "<< 500.0*(new_CC-old_CC)<<" " <<(new_clash-old_clash)<<" "<<(new_Ehbond-old_Ehbond)<<" "<<(new_Ebondlen-old_Ebondlen)<<endl;
                        dE = new_dE - old_dE;
                    //    cout<<"dE: "<<dE<<endl;
                    //    cout<<"dE_CC: "<<50.0*(new_CC - old_CC)<<endl;
        //                dE = 1.0/new_dE - 1.0/old_dE;
                    //    if(new_dE < old_dE)
                        if(dE<0.0)
                        {
                            int tmp_tmx=0;
                            for(int j=0;j<pnum;j++)
                            {
                                pointsBx[3*j+0] = tmp_mat[3*j+0];
                                pointsBx[3*j+1] = tmp_mat[3*j+1];
                                pointsBx[3*j+2] = tmp_mat[3*j+2];
                            //    bb[j]=bby[tmp_tmx];
                                tmp_tmx = tmp_tmx +1;
                            }
                            if(new_dE<best_E)
                            {
                                best_E = new_dE;
                                vector<poseCoord>().swap(best_model);
                                best_model = tmp_mat;
                            }                            
                            for(int j=0;j<pnum;j++)
                            {
                                bb[j].indn = 3*j+0;
                                bb[j].indca = 3*j+1;
                                bb[j].indc = 3*j+2; 
                                string resng = pointsrenm[j];
                                char resny = resng[0];
                                bb[j].resid = aminoid(resny);                                    
                            }     
                        //    calcsse2(bb, pnum, tmp_mat); 
                            calcssennhoc(bb,pnum,tmp_mat);
                        //    rhoC0 =rhoC01;
                    //        cout<<"KT E: "<< KT <<" "<<new_dE<<endl;
                        //    cout<<"new,old: "<< new_dE<<" "<<old_dE<<endl;
                            Aprate0 = Aprate0 +1;                             
                        //    cout<<" new_dE: "<<new_dE<<endl; 
                        }                
                        else
                        {
                        //    float tmpx=rand()/double(RAND_MAX);
                            float tmpx = randf0and1();
                            float mc_v = exp(-float(dE)/float(KT));
                            if(tmpx < mc_v)
                            {
                        //        cout<<"tmpx,mc_v,dE: "<<tmpx<<" "<<mc_v<<" "<<dE<<endl;
                        //        cout<<"new,old: "<< new_dE<<" "<<old_dE<<endl;
                        //        cout<<"KT E: "<< KT <<" "<<new_dE<<endl;
                                int tmp_tmx=0;
                                for(int j=0;j<pnum;j++)
                                {
                                    pointsBx[3*j+0] = tmp_mat[3*j+0];
                                    pointsBx[3*j+1] = tmp_mat[3*j+1];
                                    pointsBx[3*j+2] = tmp_mat[3*j+2];
                                //    bb[j]=bby[tmp_tmx];
                                    tmp_tmx = tmp_tmx +1;
                                }
                                for(int j=0;j<pnum;j++)
                                {
                                    bb[j].indn = 3*j+0;
                                    bb[j].indca = 3*j+1;
                                    bb[j].indc = 3*j+2;  
                                    string resng = pointsrenm[j];
                                    char resny = resng[0];
                                    bb[j].resid = aminoid(resny);                                       
                                }     
                            //    calcsse2(bb, pnum, tmp_mat); 
                                calcssennhoc(bb,pnum,tmp_mat);                        
                            //    rhoC0 =rhoC01;
                                Aprate0 = Aprate0 +1;                                
                            //    cout<<" XXX new_dE: "<<new_dE<<endl; 
                            } 

                        }  
        //                tmp_mat.clear(); 
            //    }

            } // the recile   
            Acprate[jjb] = float(Aprate0)/float(allnum);                                            
        }
        cout<<"ACP rate: ";
        for(int j=0;j<Acprate.size();j++)
        {
            cout<<Acprate[j]<<" ";
            if(j%15==0)
            {
                cout<<endl;
            }
        } 
    finish=clock();
    totaltime=(double)(finish-start)/CLOCKS_PER_SEC;  
    if(double(totaltime)>171000.0) break;      

}

for(int j=0;j<pnum;j++)
{
    pointsBx[3*j+0] = best_model[3*j+0];
    pointsBx[3*j+1] = best_model[3*j+1];
    pointsBx[3*j+2] = best_model[3*j+2];
//    bb[j]=bby[tmp_tmx];
} 

/*    vector<int > beg_end;
    getLinLen(randop,pointsBx,bb,beg_end);
    cout<<"beg end: "<<beg_end[0]<<" "<<beg_end[1]<<endl;

    vector<poseCoord> fra_tmp;
    for(int tt=beg_end[0];tt<=beg_end[1];tt++)
    {
        fra_tmp.push_back(pointsBx[3*tt+0]);
        fra_tmp.push_back(pointsBx[3*tt+1]);
        fra_tmp.push_back(pointsBx[3*tt+2]);
    }
    connectgap(fra_tmp);
    int tmpee=0;
    for(int tt=beg_end[0];tt<=beg_end[1];tt++)
    {
        pointsBx[3*tt+0]=fra_tmp[3*tmpee+0];
        pointsBx[3*tt+1]=fra_tmp[3*tmpee+1];
        pointsBx[3*tt+2]=fra_tmp[3*tmpee+2];
        tmpee = tmpee +1;
    }   */ 
  
    for(int i=0;i<posex.chains.size();i++)
    {
//      cout<< "chainx size: "<<posex.chains.size()<<endl;
//      cout<< "chainy size: "<<posey.chains.size()<<endl;
        Chain Chanx = posex.chains[i];
        Chain Chany;
        Chany.chaid = Chanx.chaid;
    //    fin_pose.chains.push_back(Chanx);
//      Chain Chany = posey.chains[i];
        int tmp_nx=0;
        for(int j=0;j<pnum;j++) 
        {
//          cout<<"j: "<<j<<endl;
//          cout<< "residuex size: "<< Chanx.residues.size()<<endl;
//          cout<< "residuey size: "<< Chany.residues.size()<<endl;
            Residue Resdx = Chanx.residues[j];
            Residue Resdy;
            Resdy.resseq = Resdx.resseq;
            Resdy.icode = Resdx.icode;
            Resdy.resname = Resdx.resname;
        //    fin_pose.chains[i].residues.push_back(Resdx);
//          Residue Resdy = Chany.residues[j];
//          Rot rotx,roty;
//          cout<<"Res: "<<Resdx.atoms.size()<<endl;
//          cout<<"res: "<<points_pre3[tmp_nx].res_atm.size()<<endl;
            int tt=0;
            for(int t=0;t<Resdx.atoms.size();t++)
            {
//              cout<<"t: "<<t<<endl;
//              cout<< "atomx size: "<<Resdx.atoms.size()<<endl;
//              cout<< "atomy size: "<<Resdy.atoms.size()<<endl;
                Atom Atmx=Resdx.atoms[t];
                Atom Atmy;
            //    fin_pose.chains[i].residues[j].atoms.push_back(Atmx);
//              Atom Atmy=Resdy.heavyatoms[t];                                 
                if(Atmx.atona == " CA ") 
                {
        //            fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsBx[j*3+1].x_;
                    Atmy = Atmx;
                    Atmy.xyzVect = pointsBx[j*3+1].x_;
                    Resdy.heavyatoms.push_back(Atmy);
                    Resdy.atoms.push_back(Atmy);
            //        if(j==63) cout<<"xx: "<<pointsBx[j*3+1].x_[0]<<" "<<pointsBx[j*3+1].x_[1]<<" "<<pointsBx[j*3+1].x_[2]<<endl;
            //        fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsBx[j].x_;

                } 
                if(Atmx.atona == " N  ") 
                {
    //                fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsBx[3*j].x_;
                    Atmy = Atmx;
                    Atmy.xyzVect = pointsBx[3*j].x_;
                    Resdy.heavyatoms.push_back(Atmy);
                    Resdy.atoms.push_back(Atmy);
//                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] + pointsBx[j].x_[0] - pointsBx0[j].x_[0];
//                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] + pointsBx[j].x_[1] - pointsBx0[j].x_[1];
//                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] + pointsBx[j].x_[2] - pointsBx0[j].x_[2];
                }
                if(Atmx.atona == " C  ") 
                {
//                    fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsBx[3*j+2].x_;
                    Atmy = Atmx;
                    Atmy.xyzVect = pointsBx[3*j+2].x_;
                    Resdy.heavyatoms.push_back(Atmy);
                    Resdy.atoms.push_back(Atmy);
//                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] + pointsBx[j].x_[0] - pointsBx0[j].x_[0];
//                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] + pointsBx[j].x_[1] - pointsBx0[j].x_[1];
//                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] + pointsBx[j].x_[2] - pointsBx0[j].x_[2];
                }
        /*        if(Atmx.atona == " O  ") 
                {
    //              fin_pose.chains[i].residues[j].atoms[t].xyzVect = CBp[j].Ocd;
                    float C_O_length = 1.229;
                    float CA_C_O_angle = 120.09;
                    float N_CA_C_O_diangle = 179.67;
                    vector<float> O_coord;
                    O_coord=calculateCoordinates(pointsBx[3*j].x_, pointsBx[j*3+1].x_, pointsBx[3*j+2].x_, C_O_length, CA_C_O_angle, N_CA_C_O_diangle);
                    Atmy = Atmx;
                    Atmy.xyzVect = O_coord;
                    Resdy.heavyatoms.push_back(Atmy);
                    Resdy.atoms.push_back(Atmy);
    //                fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] + pointsBx[3*j+2].x_[0] - pointsBx0[3*j+2].x_[0];
    //                fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] + pointsBx[3*j+2].x_[1] - pointsBx0[3*j+2].x_[1];
    //                fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] + pointsBx[3*j+2].x_[2] - pointsBx0[3*j+2].x_[2];
    //                fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] + pointsBx[j].x_[0] - pointsBx0[j].x_[0];
    //                fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] + pointsBx[j].x_[1] - pointsBx0[j].x_[1];
    //                fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] + pointsBx[j].x_[2] - pointsBx0[j].x_[2];                
                }          */                                                     
            }   
            tmp_nx = tmp_nx + 1; 
            Chany.residues.push_back(Resdy);
//          cout<< "tmP: "<<tmp_nx<<endl;
        }
        fin_pose.chains.push_back(Chany);
    }

    cout<< "XXXXXXXXXX"<<endl;
    return fin_pose;    
}



int main(int argc, char* argv[])
{
    char inputPDB1[600]; 
    string inputMRC;
//  poseCoords pose1,pose2;
    Model pose;

    clock_t start,finish;
    double totaltime;
    start=clock();  

//  readPDB
    strcpy(inputPDB1,argv[1]);
    readpdbstructurex(inputPDB1,pose);

    poseCoords posey;
//    readPDBcoords(inputPDB1,posey);
    cout<< "Read PDB file : "<<inputPDB1<<endl;

/*    for(int i=0;i<5;i++)
    {
        for(int j=0;j<5;j++)
        {
            cout<<randIntCustom(0,20)<<" "<<randf0and1()<<endl;
        }
    } */


    float MRC_reso;
    float mapsampling = 0.0;
    inputMRC = argv[2];
    MRC_reso = atof(argv[3]);
    mapsampling = atof(argv[4]);
//    int dm_beg= atoi(argv[5]);
//    int dm_end = atoi(argv[6]);
    char bindir[600];
    strcpy(bindir,argv[5]);
    string dist_name;
    dist_name = argv[6];
    char bindiry[600];
    strcpy(bindiry,argv[7]);

//  strcpy(inputPDB2,argv[2]);
//  readPDBcoordss(inputPDB1,pose1);
//  cout<<"000  "<<endl;
//  readPDBcoordss(inputPDB2,pose2);
//  readpdbstructurex(inputPDB2,pose2);
//  cout<< "11 "<<endl;
//  cout<<" 11 "<<endl;
    Model fin_p,fin_px,fin_py;  
    cout<<" 11 "<<endl;
    fin_p=REMC_sampley(pose,inputMRC,MRC_reso,mapsampling,dist_name,bindir,bindiry);
/*    string outnamex="11111.pdb";
    writePDBStructure(inputPDB1,fin_p,outnamex);
    cout<<"0000"<<endl;
    char cmd[1000];
    sprintf(cmd, "./Scwrl4 -i 11111.pdb -o 22222.pdb");
    system(cmd);

    cout<<"2222"<<endl;

    char rmcmd[100];
    sprintf(rmcmd, "rm -f 11111.pdb");
    system(rmcmd);

    cout<<"3333"<<endl;

    string inputPDB2="22222.pdb";
    Model poseyy;
    readpdbstructurex(inputPDB2.data(),poseyy);   
    sprintf(rmcmd, "rm -f 22222.pdb");
    system(rmcmd); */

//    string outname;
    string outnamex;
    vector<string> tmp_str,tmp_strx;
    tmp_str = string_splitx(argv[1],'.');
    tmp_strx = string_splitx(tmp_str[0],'/'); 
    string tmp_rx;
    if(tmp_strx.size()>0)
    {
        outnamex = tmp_strx[tmp_strx.size()-1] + "_11.pdb";
        tmp_rx = tmp_strx[tmp_strx.size()-1];
    }
    else
    {
        outnamex = tmp_strx[0] + "_11.pdb";
        tmp_rx = tmp_strx[0];
    }

//    string outnamex="11111.pdb";
    writePDBStructure(inputPDB1,fin_p,outnamex);
    cout<<"0000"<<endl;
//    cout<<"bindir: "<<bindir<<endl;
//    cout<<"outnamex: "<<outnamex<<endl;
//    cout<<"tmp_rx: "<<tmp_rx<<endl;
    char cmd[1000];
    sprintf(cmd, "%s/pulchra -c %s",bindir,outnamex.c_str());
    system(cmd);

    sleep(10);
    cout<<"2222"<<endl;

    char rmcmd[100];
    sprintf(rmcmd, "rm -f %s",outnamex.c_str());
    system(rmcmd);

    cout<<"3333"<<endl;

    string inputPDB2=tmp_rx+ "_11.rebuilt.pdb";
    Model poseyy;
    readpdbstructurex(inputPDB2.data(),poseyy);   
    sprintf(rmcmd, "rm -f %s",inputPDB2.c_str());
    system(rmcmd); 
    


    cout<< "444"<<endl;
    fin_py = poseyy;  
//  fin_px=REMC_samplex(fin_p,pose2);
//  fin_py=REMC_samplexy(fin_px,pose2);
//  writePDBcoordss(inputPDB1, fin_p);
    string outname;
/*    vector<string> tmp_str,tmp_strx;
    tmp_str = string_splitx(argv[1],'.');
    tmp_strx = string_splitx(tmp_str[0],'/'); */

    if(tmp_strx.size()>0)
    {
        outname = tmp_strx[tmp_strx.size()-1] + "_outp00x03.pdb";
    }
    else
    {
        outname = tmp_strx[0] + "_outp00x03.pdb";
    }
    
    cout<<"out: "<<outname<<endl;
    writePDBStructure(inputPDB1,fin_py,outname);
    finish=clock();
    totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
    cout<<"\nrunning time: "<<totaltime<<" s！"<<endl;  
    return 0;
}