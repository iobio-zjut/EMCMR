#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <math.h>
#include <time.h>

using namespace std;

struct press
{
	string ss_type;
	string res_name;
	int  ss_ord;
};

void readssfile(string filename, vector<press> &presstype)
{
    ifstream ReadFile;
    ReadFile.open(filename.c_str());
    if(!ReadFile.is_open())
    {
        cout<<"error to open seq.dat (ss) the file"<<endl;
        return ;
    }   
    string atomline; 
    while(!ReadFile.eof())
    {
        getline(ReadFile,atomline);
        if(atomline.length()>0)
        {
	        press pret;
	        pret.ss_ord = (stoi(atomline.substr(0,6)));
	        pret.res_name = (atomline.substr(6,7));
//	        cout<<"1"<<atomline.substr(6,7)<<2<<endl;
	        pret.ss_type = (atomline.substr(13,4));
	        presstype.push_back(pret);
        }
    }
}
int main(int argc, char* argv[])
{
	string filename;
	string bindir = argv[2];
	string seq_name = argv[1];
	filename = bindir + "/" + seq_name;
	cout<<filename<<endl;
	vector<press> presstype;
	readssfile(filename,presstype);
	for(int i=0;i<presstype.size();i++)
	{
		cout<<"ss_ord: "<<presstype[i].ss_ord<<" "<<presstype[i].res_name<<" "<<presstype[i].ss_type<<endl;
	}
}