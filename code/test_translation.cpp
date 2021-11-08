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
float Distance_pointy(vector<vector<float> > A,vector<vector<float> > B)
{

}
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

Model REMC_sampley(Model posex, string inputM,float MRC_R,float mapsp)
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
    theDensityMap.pCA_MASK =7.0;
    theDensityMap.pforce_apix_on_map_load_= 0.0;
    theDensityMap.pWINDOW_ =1;
    theDensityMap.pscore_window_context_ = false ;
    theDensityMap.premap_symm_ = false ;
    theDensityMap.pde_edensity = false;
    theDensityMap.pnkbins_ = 0 ;   

    theDensityMap.readMRCandResize(inputM,MRC_R,mapsp);
    theDensityMap.writeMRC("xx1.mrc");
//  theDensityMap.calcRhoC(posey,4.0,rhoC,rhoMask);
//    theDensityMap.calcRhoCx(posex,4.0,rhoC,rhoMask);
//  theDensityMap.calcRhoCy(pose,4.0,rhoC,rhoMask);
//    float CC;
//  CC=theDensityMap.getRSCC(rhoC,rhoMask);
 //   CC=theDensityMap.getRSCCX(rhoC,rhoMask);
 //   cout<<"CC: "<<CC<<endl; 
    cout<<"origin: "<<theDensityMap.origin[0]<<" "<<theDensityMap.origin[1]<<" "<<theDensityMap.origin[2]<<endl;;

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
    vector<poseCoord> best_model;
    float best_E=10.0;

    vector<vector<float>> atm_idx(pnum,vector<float>(3,0.0));


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
//            k=(M_PI/effReso)*(M_PI/effReso); 
//            a= sig_j.s_mass;
//            C= a*pow(k/3.1415926,1.5);                
            k = sig_j.k( theDensityMap.effectiveB );
            C = sig_j.C( k );
            if ( C < 1e-6 ) continue;   

            cartX1 = pointsBx[i].x_;
            cartX1[0] = cartX1[0] - theDensityMap.origin[0];
            cartX1[1] = cartX1[1] - theDensityMap.origin[1];
            cartX1[2] = cartX1[2] - theDensityMap.origin[2];

            MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);

            // the location of atom in grid ?
            atm_idx[tmp_i][0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] + 1) , (double)theDensityMap.grid[0]);
            atm_idx[tmp_i][1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] + 1) , (double)theDensityMap.grid[1]);
            atm_idx[tmp_i][2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] + 1) , (double)theDensityMap.grid[2]);   
    //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
    //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
    //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
    //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
            for(int z=1;z<=theDensityMap.density.u3();z++)
            {
                atm_jx[2] = z;
                del_ijx[2] =(atm_idx[tmp_i][2]-atm_jx[2])/theDensityMap.grid[2];
                if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                del_ijx[0] = del_ijx[1] = 0.0;
                vector<float> frac_tmpz;
                MatrixTimesTransVector(theDensityMap.f2c,del_ijx,frac_tmpz);
                if(square_len(frac_tmpz)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                    
                for(int y=1;y<=theDensityMap.density.u2();y++)
                {
                    atm_jx[1] = y;
                    del_ijx[1] = (atm_idx[tmp_i][1] - atm_jx[1])/theDensityMap.grid[1];
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
                        del_ijx[0] = (atm_idx[tmp_i][0] - atm_jx[0])/theDensityMap.grid[0];
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
        if(rhoC0[x]>max_rhc) max_rhc = rhoC0[x] ;
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
    for(int x=1;x<=theDensityMap.density.u1();++x)
    {
       for(int y=1;y<=theDensityMap.density.u2();++y)
        {
            for(int z=1;z<=theDensityMap.density.u3();++z)
            {
                clc_x2 = rhoC0(x,y,z);
//                clc_x2 = tmp_rhc;
                obs_x2 = theDensityMap.density(x,y,z);
//                obs_x2 = tmp_den;
                eps_x2 = 1.0-inv_rho_mask0(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                sumCO_i2 += eps_x2*clc_x2*obs_x2;
                sumO_i2  += eps_x2*obs_x2;
                sumO2_i2 += eps_x2*obs_x2*obs_x2;
                sumC_i2  += eps_x2*clc_x2;
                sumC2_i2 += eps_x2*clc_x2*clc_x2;
                vol_i2   += eps_x2;
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
    varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
    varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 );
    if ( varC_i2 == 0 || varO_i2 == 0 ) {
        CC_i2 = 0;
    } else {
        CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
    }   
    best_E=1.0-CC_i2;  
    cout<<"best_E: "<<best_E<<endl;  
    best_model = pointsBx;
//    mass_pdb[0] = mass_pdb[0]/allpdbmass;
//    mass_pdb[1] = mass_pdb[1]/allpdbmass;
//    mass_pdb[2] = mass_pdb[2]/allpdbmass;
//    mass_mrc[0] = mass_mrc[0]/allmrcmass;
//    mass_mrc[1] = mass_mrc[1]/allmrcmass;
//    mass_mrc[2] = mass_mrc[2]/allmrcmass; 
//    mass_mrc[0] = theDensityMap.density.u1()/2.0 + (theDensityMap.origin[0]);//*(1.0-1.0/2.0);
//    mass_mrc[1] = theDensityMap.density.u2()/2.0 + (theDensityMap.origin[1]);//*(1.0-1.0/2.0);
//    mass_mrc[2] = theDensityMap.density.u3()/2.0 + (theDensityMap.origin[2]);//*(1.0-1.0/2.0);        
    mass_mrc[0] = ((float(theDensityMap.cellDimensions[0])/float(theDensityMap.density.u1()))*float(theDensityMap.grid[0]))/2.0 + (theDensityMap.origin[0]);//*(1.0-1.0/2.0);
    mass_mrc[1] = ((float(theDensityMap.cellDimensions[1])/float(theDensityMap.density.u2()))*float(theDensityMap.grid[1]))/2.0 + (theDensityMap.origin[1]);//*(1.0-1.0/2.0);
    mass_mrc[2] = ((float(theDensityMap.cellDimensions[2])/float(theDensityMap.density.u3()))*float(theDensityMap.grid[2]))/2.0 + (theDensityMap.origin[2]);//*(1.0-1.0/2.0);    
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
//    mass_mrc[0] = mass_mrc[0]/theDensityMap.grid[0];
//    mass_mrc[1] = mass_mrc[1]/theDensityMap.grid[1];
//    mass_mrc[2] = mass_mrc[2]/theDensityMap.grid[2];
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
    mass_mrc0 = mass_mrc;
//    MatrixTimesTransVector(theDensityMap.f2c,mass_mrc,mass_mrc0);
    cout<<"mass_mrc: "<<mass_mrc0[0]<<" "<<mass_mrc0[1]<<" "<<mass_mrc0[2]<<endl;
//    vector<float> f_origin(3,0.0);
//    f_origin[0] = theDensityMap.origin[0];
//    f_origin[1] = theDensityMap.origin[1];
//    f_origin[2] = theDensityMap.origin[2];
//    vector<float> p_origin(3,0.0);
//    MatrixTimesTransVector(theDensityMap.f2c,f_origin,p_origin);
//    mass_mrc0[0] = mass_mrc0[0] - p_origin[0];
//    mass_mrc0[1] = mass_mrc0[1] - p_origin[1];
//    mass_mrc0[2] = mass_mrc0[2] - p_origin[2];
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
    vector<poseCoord> pointsBxj;
    pointsBxj= pointsBx;
    for(int j=0;j<pointsBx.size();j++) // not last N
    {
        pointsBxj[j].x_[0] = tran0[0] + pointsBxj[j].x_[0];
        pointsBxj[j].x_[1] = tran0[1] + pointsBxj[j].x_[1];
        pointsBxj[j].x_[2] = tran0[2] + pointsBxj[j].x_[2];
//        tmp3=tmp3+1;
    }  
//    pointsBx = pointsBxj;
    clock_t start,finish;
    double totaltime;
    start=clock();  
    rhoC0.dimension(theDensityMap.density.u1() , theDensityMap.density.u2() , theDensityMap.density.u3());
//    rhoC01.dimension(theDensityMap.density.u1() , theDensityMap.density.u2() , theDensityMap.density.u3());
    inv_rho_mask0.dimension(theDensityMap.density.u1() , theDensityMap.density.u2() , theDensityMap.density.u3());
//    inv_rho_mask1.dimension(theDensityMap.density.u1() , theDensityMap.density.u2() , theDensityMap.density.u3());
    for ( int t=0; t<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++t ) {
        rhoC0[t]=0.0;
        inv_rho_mask0[t]=1.0;
    }
    del_ijx=vector<float>(3,0.0);
    atm_jx=vector<float>(3,0.0); 
//    string elt_i; 
    tmp_i=0;
    bond0max=-1000.0;bond1max=-1000.0;bond2max=-1000.0;
    bond0min=1000.0;bond1min=1000.0;bond2min=1000.0;  
//    cout<<"FF"<<endl;      
    for(int i=0; i<pointsBxj.size();i++)
    {
//        cout<<"YYYY: "<<pointsBx[i].elt_<<endl;;
//        if(pointsBx[i].elt_ == "CA")
//        {
         if(pointsBxj[i].elt_ == "CA")
        {           
//            cout<<"RRRR: "<<i;
            vector<float> cartX1;
            vector<float> fracX1;
            elt_i = pointsBxj[i].elt_;
            elt_i = elt_i[0];
            OneGaussianScattering sig_j = get_A( elt_i );
//            k=(M_PI/effReso)*(M_PI/effReso); 
//            a= sig_j.s_mass;
//            C= a*pow(k/3.1415926,1.5);            
            k = sig_j.k( theDensityMap.effectiveB );
            C = sig_j.C( k );
            if ( C < 1e-6 ) continue;   

            cartX1 = pointsBxj[i].x_;  
            cartX1[0] = cartX1[0] - theDensityMap.origin[0];
            cartX1[1] = cartX1[1] - theDensityMap.origin[1];
            cartX1[2] = cartX1[2] - theDensityMap.origin[2];                     
            MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);

            // the location of atom in grid ?
            vector<float> atm_idxt(3,0.0);
            atm_idxt[0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0]  + 1) , (double)theDensityMap.grid[0]);
            atm_idxt[1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1]  + 1) , (double)theDensityMap.grid[1]);
            atm_idxt[2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2]  + 1) , (double)theDensityMap.grid[2]);   
    //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
    //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
    //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
    //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
            for(int z=1;z<=theDensityMap.density.u3();z++)
            {
                atm_jx[2] = z;
                del_ijx[2] =(atm_idxt[2]-atm_jx[2])/theDensityMap.grid[2];
                if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                del_ijx[0] = del_ijx[1] = 0.0;
                vector<float> frac_tmpz;
                MatrixTimesTransVector(theDensityMap.f2c,del_ijx,frac_tmpz);
                if(square_len(frac_tmpz)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                    
                for(int y=1;y<=theDensityMap.density.u2();y++)
                {
                    atm_jx[1] = y;
                    del_ijx[1] = (atm_idxt[1] - atm_jx[1])/theDensityMap.grid[1];
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
                        del_ijx[0] = (atm_idxt[0] - atm_jx[0])/theDensityMap.grid[0];
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
/*    float max_den = 0.0;
    float max_rhc = 0.0;
    for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
//        cout<<theDensityMap.density[x]<<" ";
        if(theDensityMap.density[x]>max_den) max_den = theDensityMap.density[x] ;
        if(rhoC0[x]>max_rhc) max_rhc =rhoC0[x] ;
    } */
    mass_pdb=vector<float>(3,0.0);
    mass_mrc=vector<float>(3,0.0);
    coor_pdb=vector<float>(3,0.0);  
    allpdbmass=0.0;
    allmrcmass=0.0;
    sumC_i2=0.0, sumO_i2=0.0, sumCO_i2=0.0, vol_i2=0.0, CC_i2=0.0;
    sumO2_i2=0.0, sumC2_i2=0.0, varC_i2=0.0, varO_i2=0.0;
    clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;    
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
                sumCO_i2 += eps_x2*clc_x2*obs_x2;
                sumO_i2  += eps_x2*obs_x2;
                sumO2_i2 += eps_x2*obs_x2*obs_x2;
                sumC_i2  += eps_x2*clc_x2;
                sumC2_i2 += eps_x2*clc_x2*clc_x2;
                vol_i2   += eps_x2;
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
    varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
    varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 );
    if ( varC_i2 == 0 || varO_i2 == 0 ) {
        CC_i2 = 0;
    } else {
        CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
    }   
    
    finish=clock();
    totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
    cout<<"\nrunning time: "<<totaltime<<" s！"<<endl; 


    float E_2=1.0-CC_i2; 
    cout<<"E_2: "<<E_2<<endl;       
    if(E_2<=best_E)
    {
        vector<poseCoord>().swap(pointsBx);
        pointsBx=pointsBxj;
    }
    cout<<"CCCCC"<<endl;    
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
    float mapspx=4.0;
    theDensityMapx.readMRCandResize(inputM,MRC_R,mapspx); 
    theDensityMapx.writeMRC_test("xx5.mrc"); 
    cout<<"density,x y z: "<<theDensityMapx.density.u1()<<" "<<theDensityMapx.density.u2()<<" "<<theDensityMapx.density.u3()<<endl;
    float new_CC1=0.0;
    float old_CC1=0.0;
    float new_CC2=0.0;     
//    vector<poseCoord> best_model;
//    float best_E=10.0;
    ObjexxFCL::FArray3D< float > rhoC0x;
    ObjexxFCL::FArray3D< float > inv_rho_mask0x;
    rhoC0x.dimension(theDensityMapx.density.u1() , theDensityMapx.density.u2() , theDensityMapx.density.u3());
    inv_rho_mask0x.dimension(theDensityMapx.density.u1() , theDensityMapx.density.u2() , theDensityMapx.density.u3());
    // find initial position for A and B
    float init0max=-10000.0,init1max=-10000.0,init2max=-10000.0;
    float init0min=10000.0,init1min=10000.0,init2min=10000.0;       
//    best_E=10.0;
    coor_pdb=vector<float>(3,0.0);   
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
            k=(M_PI/effReso)*(M_PI/effReso); 
            a= sig_j.s_mass;
            C= a*pow(k/3.1415926,1.5);

//            k = sig_j.k( theDensityMapx.effectiveB );
//            C = sig_j.C( k );
//            if ( C < 1e-6 ) continue;   

            cartX1 = pointsBx[i].x_;   
            cartX1[0] = cartX1[0] - theDensityMapx.origin[0];
            cartX1[1] = cartX1[1] - theDensityMapx.origin[1];
            cartX1[2] = cartX1[2] - theDensityMapx.origin[2];                      
            MatrixTimesTransVector(theDensityMapx.c2f,cartX1,fracX1);

            // the location of atom in grid ?
            atm_idx[tmp_i][0] = pos_mod (double(fracX1[0]*theDensityMapx.grid[0]  + 1) , (double)theDensityMapx.grid[0]);
            atm_idx[tmp_i][1] = pos_mod (double(fracX1[1]*theDensityMapx.grid[1]  + 1) , (double)theDensityMapx.grid[1]);
            atm_idx[tmp_i][2] = pos_mod (double(fracX1[2]*theDensityMapx.grid[2]  + 1) , (double)theDensityMapx.grid[2]);   
    //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
    //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
    //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
    //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
            for(int z=1;z<=theDensityMapx.density.u3();z++)
            {
                atm_jx[2] = z;
                del_ijx[2] =(atm_idx[tmp_i][2]-atm_jx[2])/theDensityMapx.grid[2];
                if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                del_ijx[0] = del_ijx[1] = 0.0;
                vector<float> frac_tmpz;
                MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,frac_tmpz);
                if(square_len(frac_tmpz)> (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;                    
                for(int y=1;y<=theDensityMapx.density.u2();y++)
                {
                    atm_jx[1] = y;
                    del_ijx[1] = (atm_idx[tmp_i][1] - atm_jx[1])/theDensityMapx.grid[1];
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
                        del_ijx[0] = (atm_idx[tmp_i][0] - atm_jx[0])/theDensityMapx.grid[0];
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
    for(int x=1;x<=theDensityMapx.density.u1();x++)
    {
       for(int y=1;y<=theDensityMapx.density.u2();y++)
        {
            for(int z=1;z<=theDensityMapx.density.u3();z++)
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
    cout<<"new_CC2: "<<new_CC2<<endl;
    best_E = new_CC2;
    vector<poseCoord>().swap(best_model);
    best_model = pointsBx;    
//    cout<<"DDDD"<<endl;

    cout<<"rotate 0 and 180"<<endl;
    // rotate 60 and 180
//    vector<poseCoord>().swap(pointsBx);
//    pointsBx = best_model;  
    int tmp3=0;
    float init_E0=0.0;
    float init_E1=0.0;  
    vector<poseCoord> init_mat0;
    vector<poseCoord> init_mat1;    
    for(int kkk=0;kkk<6;kkk++) // 2
    {
        vector<poseCoord > point_mat;  
        float zang=180.0;
        float zang_rd= zang*M_PI/180.0;            
        if(kkk%3==0) // rotate around Z axis
        { 
            vector<poseCoord>().swap(point_mat);
            coor_pdb=vector<float>(3,0.0);
            for(int j=0;j<pnum;j++)
            {
                point_mat.push_back(pointsBx[3*j+0]);
                point_mat.push_back(pointsBx[3*j+1]);
                point_mat.push_back(pointsBx[3*j+2]);
                coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
                coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
                coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];                     
            }              
            coor_pdb[0] = coor_pdb[0]/(pnum);
            coor_pdb[1] = coor_pdb[1]/(pnum);
            coor_pdb[2] = coor_pdb[2]/(pnum);             
            float asin_theta=2.0*randf0and1()-1.0;
            float acos_theta=sqrt(1.0-asin_theta*asin_theta);
            float apha=2.0*PI*randf0and1();
            float awx=0.0;
            float awy=0.0;
            float awz=1.0;
            // Translation Vector
            float t0=0.0;
            float angle_rotategg=zang_rd; // rotate angle  
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
            axyz[0]= coor_pdb[0];  
            axyz[1]= coor_pdb[1];
            axyz[2]= coor_pdb[2];
                
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

/*            for(int j=0;j<pnum;j++)
            {
                point_mat.push_back(best_model[3*j+0]);
                point_mat.push_back(best_model[3*j+1]);
                point_mat.push_back(best_model[3*j+2]);
            }  
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
            vector<poseCoord> point_matz;
            point_matz = point_mat;
            for(int j=0;j<pnum;j++)
            {
                point_mat[3*j+0].x_[0]= point_matz[3*j+0].x_[0]*cos(zang_rd)+point_matz[3*j+0].x_[1]*(-sin(zang_rd));
                point_mat[3*j+0].x_[1]= point_matz[3*j+0].x_[0]*sin(zang_rd)+point_matz[3*j+0].x_[1]*(cos(zang_rd));
                point_mat[3*j+0].x_[2]= point_matz[3*j+0].x_[2];  
                point_mat[3*j+1].x_[0]= point_matz[3*j+1].x_[0]*cos(zang_rd)+point_matz[3*j+1].x_[1]*(-sin(zang_rd));
                point_mat[3*j+1].x_[1]= point_matz[3*j+1].x_[0]*sin(zang_rd)+point_matz[3*j+1].x_[1]*(cos(zang_rd));
                point_mat[3*j+1].x_[2]= point_matz[3*j+1].x_[2]; 
                point_mat[3*j+2].x_[0]= point_matz[3*j+2].x_[0]*cos(zang_rd)+point_matz[3*j+2].x_[1]*(-sin(zang_rd));
                point_mat[3*j+2].x_[1]= point_matz[3*j+2].x_[0]*sin(zang_rd)+point_matz[3*j+2].x_[1]*(cos(zang_rd));
                point_mat[3*j+2].x_[2]= point_matz[3*j+2].x_[2];               
            }    */        
        }
        if(kkk%3==1)  //rotate around X axis
        {
            vector<poseCoord>().swap(point_mat);
            coor_pdb=vector<float>(3,0.0);
            for(int j=0;j<pnum;j++)
            {
                point_mat.push_back(pointsBx[3*j+0]);
                point_mat.push_back(pointsBx[3*j+1]);
                point_mat.push_back(pointsBx[3*j+2]);
                coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
                coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
                coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];                     
            }              
            coor_pdb[0] = coor_pdb[0]/(pnum);
            coor_pdb[1] = coor_pdb[1]/(pnum);
            coor_pdb[2] = coor_pdb[2]/(pnum);             
            float asin_theta=2.0*randf0and1()-1.0;
            float acos_theta=sqrt(1.0-asin_theta*asin_theta);
            float apha=2.0*PI*randf0and1();
            float awx=1.0;
            float awy=0.0;
            float awz=0.0;
            // Translation Vector
            float t0=0.0;
            float angle_rotategg=zang_rd; // rotate angle  
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
            axyz[0]= coor_pdb[0];  
            axyz[1]= coor_pdb[1];
            axyz[2]= coor_pdb[2];
                
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
/*            vector<poseCoord>().swap(point_mat);
            for(int j=0;j<pnum;j++)
            {
                point_mat.push_back(best_model[3*j+0]);
                point_mat.push_back(best_model[3*j+1]);
                point_mat.push_back(best_model[3*j+2]);
            }              
            vector<poseCoord> point_matx;
            point_matx = point_mat;    
            for(int j=0;j<pnum;j++)
            {
                point_mat[3*j+0].x_[0]= point_matx[3*j+0].x_[0];
                point_mat[3*j+0].x_[1]= point_matx[3*j+0].x_[1]*cos(zang_rd)+point_matx[3*j+0].x_[2]*(-sin(zang_rd));
                point_mat[3*j+0].x_[2]= point_matx[3*j+0].x_[1]*sin(zang_rd)+point_matx[3*j+0].x_[2]*cos(zang_rd);  
                point_mat[3*j+1].x_[0]= point_matx[3*j+1].x_[0];
                point_mat[3*j+1].x_[1]= point_matx[3*j+1].x_[1]*cos(zang_rd)+point_matx[3*j+1].x_[2]*(-sin(zang_rd));
                point_mat[3*j+1].x_[2]= point_matx[3*j+1].x_[1]*sin(zang_rd)+point_matx[3*j+1].x_[2]*cos(zang_rd); 
                point_mat[3*j+2].x_[0]= point_matx[3*j+2].x_[0];
                point_mat[3*j+2].x_[1]= point_matx[3*j+2].x_[1]*cos(zang_rd)+point_matx[3*j+2].x_[2]*(-sin(zang_rd));
                point_mat[3*j+2].x_[2]= point_matx[3*j+2].x_[1]*sin(zang_rd)+point_matx[3*j+2].x_[2]*cos(zang_rd);               
            }  */
        } 
        if(kkk%3==2)  //rotate around Y axis
        {
            vector<poseCoord>().swap(point_mat);
            coor_pdb=vector<float>(3,0.0);
            for(int j=0;j<pnum;j++)
            {
                point_mat.push_back(pointsBx[3*j+0]);
                point_mat.push_back(pointsBx[3*j+1]);
                point_mat.push_back(pointsBx[3*j+2]);
                coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
                coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
                coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];                     
            }              
            coor_pdb[0] = coor_pdb[0]/(pnum);
            coor_pdb[1] = coor_pdb[1]/(pnum);
            coor_pdb[2] = coor_pdb[2]/(pnum);             
            float asin_theta=2.0*randf0and1()-1.0;
            float acos_theta=sqrt(1.0-asin_theta*asin_theta);
            float apha=2.0*PI*randf0and1();
            float awx=0.0;
            float awy=1.0;
            float awz=0.0;
            // Translation Vector
            float t0=0.0;
            float angle_rotategg=zang_rd; // rotate angle  
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
            axyz[0]= coor_pdb[0];  
            axyz[1]= coor_pdb[1];
            axyz[2]= coor_pdb[2];
                
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
/*            vector<poseCoord>().swap(point_mat);
            for(int j=0;j<pnum;j++)
            {
                point_mat.push_back(best_model[3*j+0]);
                point_mat.push_back(best_model[3*j+1]);
                point_mat.push_back(best_model[3*j+2]);
            }             
            vector<poseCoord> point_maty;
            point_maty = point_mat;    
            for(int j=0;j<pnum;j++)
            {
                point_mat[3*j+0].x_[0]= point_maty[3*j+0].x_[0]*cos(zang_rd)+point_maty[3*j+0].x_[2]*sin(zang_rd);
                point_mat[3*j+0].x_[1]= point_maty[3*j+0].x_[1];
                point_mat[3*j+0].x_[2]= point_maty[3*j+0].x_[0]*(-sin(zang_rd))+point_maty[3*j+0].x_[2]*cos(zang_rd);  
                point_mat[3*j+1].x_[0]= point_maty[3*j+1].x_[0]*cos(zang_rd)+point_maty[3*j+1].x_[2]*sin(zang_rd);
                point_mat[3*j+1].x_[1]= point_maty[3*j+1].x_[1];
                point_mat[3*j+1].x_[2]= point_maty[3*j+1].x_[0]*(-sin(zang_rd))+point_maty[3*j+1].x_[2]*cos(zang_rd); 
                point_mat[3*j+2].x_[0]= point_maty[3*j+2].x_[0]*cos(zang_rd)+point_maty[3*j+2].x_[2]*sin(zang_rd);
                point_mat[3*j+2].x_[1]= point_maty[3*j+2].x_[1];
                point_mat[3*j+2].x_[2]= point_maty[3*j+2].x_[0]*(-sin(zang_rd))+point_maty[3*j+2].x_[2]*cos(zang_rd);               
            }  */
        }        
//        best_E=10.0;
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
                k=(M_PI/effReso)*(M_PI/effReso); 
                a= sig_j.s_mass;
                C= a*pow(k/3.1415926,1.5);                
//                k = sig_j.k( theDensityMapx.effectiveB );
//                C = sig_j.C( k );
//                if ( C < 1e-6 ) continue;   

                cartX1 = point_mat[i].x_;
                cartX1[0] = cartX1[0] - theDensityMapx.origin[0];
                cartX1[1] = cartX1[1] - theDensityMapx.origin[1];
                cartX1[2] = cartX1[2] - theDensityMapx.origin[2];                             
                MatrixTimesTransVector(theDensityMapx.c2f,cartX1,fracX1);

                // the location of atom in grid ?
                atm_idx[tmp_i][0] = pos_mod (double(fracX1[0]*theDensityMapx.grid[0]  + 1) , (double)theDensityMapx.grid[0]);
                atm_idx[tmp_i][1] = pos_mod (double(fracX1[1]*theDensityMapx.grid[1]  + 1) , (double)theDensityMapx.grid[1]);
                atm_idx[tmp_i][2] = pos_mod (double(fracX1[2]*theDensityMapx.grid[2]  + 1) , (double)theDensityMapx.grid[2]);   
        //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
        //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
        //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
        //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
                for(int z=1;z<=theDensityMapx.density.u3();z++)
                {
                    atm_jx[2] = z;
                    del_ijx[2] =(atm_idx[tmp_i][2]-atm_jx[2])/theDensityMapx.grid[2];
                    if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                    if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                    del_ijx[0] = del_ijx[1] = 0.0;
                    vector<float> frac_tmpz;
                    MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,frac_tmpz);
                    if(square_len(frac_tmpz)> (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;                    
                    for(int y=1;y<=theDensityMapx.density.u2();y++)
                    {
                        atm_jx[1] = y;
                        del_ijx[1] = (atm_idx[tmp_i][1] - atm_jx[1])/theDensityMapx.grid[1];
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
                            del_ijx[0] = (atm_idx[tmp_i][0] - atm_jx[0])/theDensityMapx.grid[0];
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
        for(int x=1;x<=theDensityMapx.density.u1();x++)
        {
           for(int y=1;y<=theDensityMapx.density.u2();y++)
            {
                for(int z=1;z<=theDensityMapx.density.u3();z++)
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
        if(new_CC2<=best_E) 
        {
            best_E = new_CC2;
            best_model = point_mat;
        }
        old_dE= new_CC2; 
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
        //    old_dE=new_CC2;   
        //    old_dE = old_CC1 ;// + old_clash + 0.1*d1/d2;// + olddis_p_m;
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
                    k=(M_PI/effReso)*(M_PI/effReso); 
                    a= sig_j.s_mass;
                    C= a*pow(k/3.1415926,1.5);                    
//                    k = sig_j.k( theDensityMapx.effectiveB );
//                    C = sig_j.C( k );
//                    if ( C < 1e-6 ) continue;   

                    cartX1 = fin_mat[i].x_;  
                    cartX1[0] = cartX1[0] - theDensityMapx.origin[0];
                    cartX1[1] = cartX1[1] - theDensityMapx.origin[1];
                    cartX1[2] = cartX1[2] - theDensityMapx.origin[2];                               
                    MatrixTimesTransVector(theDensityMapx.c2f,cartX1,fracX1);

                    // the location of atom in grid ?
                    atm_idx[tmp_i][0] = pos_mod (double(fracX1[0]*theDensityMapx.grid[0]  + 1) , (double)theDensityMapx.grid[0]);
                    atm_idx[tmp_i][1] = pos_mod (double(fracX1[1]*theDensityMapx.grid[1]  + 1) , (double)theDensityMapx.grid[1]);
                    atm_idx[tmp_i][2] = pos_mod (double(fracX1[2]*theDensityMapx.grid[2]  + 1) , (double)theDensityMapx.grid[2]);   
            //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
            //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
            //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
            //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
                    for(int z=1;z<=theDensityMapx.density.u3();z++)
                    {
                        atm_jx[2] = z;
                        del_ijx[2] =(atm_idx[tmp_i][2]-atm_jx[2])/theDensityMapx.grid[2];
                        if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                        if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                        del_ijx[0] = del_ijx[1] = 0.0;
                        vector<float> frac_tmpz;
                        MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,frac_tmpz);
                        if(square_len(frac_tmpz)> (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;                    
                        for(int y=1;y<=theDensityMapx.density.u2();y++)
                        {
                            atm_jx[1] = y;
                            del_ijx[1] = (atm_idx[tmp_i][1] - atm_jx[1])/theDensityMapx.grid[1];
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
                                del_ijx[0] = (atm_idx[tmp_i][0] - atm_jx[0])/theDensityMapx.grid[0];
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
            for(int x=1;x<=theDensityMapx.density.u1();x++)
            {
               for(int y=1;y<=theDensityMapx.density.u2();y++)
                {
                    for(int z=1;z<=theDensityMapx.density.u3();z++)
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
            new_dE=1.0-CC_i2; 
    //        new_dE = new_CC1 ;                                   
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
        //        new_CC2 = new_CC1;
                old_dE = new_dE;
                 
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
                    old_dE = new_dE;
                //    new_CC2 = new_CC1;
                    cout<<"old_CC1,new_CC1:XX "<<old_dE<<" "<<new_dE<<endl;                    
                }                
            }           
        }
/*        if(kkk==0)
        {
            init_E0=best_E;
            vector<poseCoord>().swap(init_mat0);
            init_mat0 = best_model;
        }
        if(kkk==1)
        {
            init_E1=best_E;
            vector<poseCoord>().swap(init_mat1);
            init_mat1 = best_model;            
        }    */    
    }
    vector<poseCoord>().swap(pointsBx);
    pointsBx=best_model;



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

//  strcpy(inputPDB2,argv[2]);
//  readPDBcoordss(inputPDB1,pose1);
//  cout<<"000  "<<endl;
//  readPDBcoordss(inputPDB2,pose2);
//  readpdbstructurex(inputPDB2,pose2);
//  cout<< "11 "<<endl;
//  cout<<" 11 "<<endl;
    Model fin_p,fin_px,fin_py;  
    cout<<" 11 "<<endl;
    fin_p=REMC_sampley(pose,inputMRC,MRC_reso,mapsampling);
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
