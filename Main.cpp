#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <numeric>

#include "Handy.h"


/*==========================================================================
This parses the bam file to generate values by transcript
Reports:
Reads
Unique reads
Pairs where both ends align 
Duplicates
Unique Duplicates


==========================================================================*/

unsigned int checkErrors();
void displayHelp();
void somethingsGoneWrong(unsigned int);
void processFastqFile();
void writeResults();

//Variables taken from the inputs
string fastqFileName="";
string outFileName="";
int nReads=0;
int seed=4503;
int merLength=5;
double cutOff=0.0001;
int lineCt=0;//how many lines were processed
int totalReads=0;
int samplingPercent=20;
bool firstMerOnly=false;


map<string, vector<int> > merMap;
map<string, int > merReadMap;
vector<int> emptyVector;
map<string, vector<int> > baseMap;
map<string, vector<double> > obsMap;
vector<int>totalMersAtPosition;

using namespace std;

/*
 *Useage: 

 * ./GetMersFastq -out ./testChris -fastq /seq/mbrd/Software/Pipeline/output/2013_11_21_GeorgiaChimera_H79AAADXX/ERP_RZ_noPNK_ATGTTGTA/Reads/left_read.fastq -mer 10
*/
int main(int argc, char* argv[]) 
{
		
	//Intialize library of useful things
	Handy h(0);
	
	int optind=1;

	while ((optind < argc) && (argv[optind][0]=='-')) 
	{	
        string sw = argv[optind];
				
		if (sw=="-h") 
		{	
            optind++;
			displayHelp();
			return 1;
        }
		
		else if(optind >=  argc-1 && sw != "-firstMer")
		{
			cerr<<"Your final parameter, "<<sw<<" is missing a value."<<endl;
			return 1;
		}

		else if (h.cmpStringNoCase(sw, "-fastq")==1)
		{	
            optind++;
			fastqFileName = argv[optind];		
			optind++;
        }
		
		else if (h.cmpStringNoCase(sw, "-out")==1)
		{	
            optind++;
			outFileName = argv[optind];		
			optind++;
        }
		
		else if (h.cmpStringNoCase(sw, "-cutOff")==1)
		{	
            optind++;
			cutOff = h.getDoubleFromString(argv[optind]);		
			optind++;
        }
		
		else if (h.cmpStringNoCase(sw, "-mer")==1)
		{	
            optind++;
			merLength = h.getIntFromString(argv[optind]);		
			optind++;
        }
		
		else if (h.cmpStringNoCase(sw, "-sampling_percent")==1)
		{	
            optind++;
			samplingPercent= h.getIntFromString(argv[optind]);		
			optind++;
        }
		
		else if (h.cmpStringNoCase(sw, "-firstMer")==1)
		{	
            optind++;
			firstMerOnly = true;		
			
        }
		
		else
		{
			cerr<<"Main: Unknown parameter:"<<sw<<endl;
			return 1;
		}
	}	
	
	checkErrors();

	srand (4503);
	
	processFastqFile();
	cout<<"Current size:"<<merMap.size();

	cout<<"Done\n";
	
	writeResults();
	cout<<"Unique Reads Found."<<endl;

		
}


void processFastqFile()
{
	Handy h(0);

	ifstream FASTQ(fastqFileName.c_str(), ios::in);	 
	string line;
	string seq;
	vector <int> emptyVector2;
	int longestLine=0;
	lineCt=0;

    // ensure that the file can be opened
	if ( !FASTQ )
	{
		cerr << "ERROR!! FASTQ file: " << fastqFileName << " could not be opened" << endl;
		exit (1);
	}
  
	cout<<"FASTQ file opened, about to begin reading in."<<endl;
	
	while (getline(FASTQ, line , '\n' )) 
	{
		if(lineCt%4==1)
		{
			if((int)line.length()>longestLine)
			{
				longestLine=(int)line.length();
			}
		}
		++lineCt;
	}
	FASTQ.close();

	cout<<"Longest line found:"<<longestLine<<endl;
	
	lineCt=0;
	
	ifstream FASTQ1(fastqFileName.c_str(), ios::in);	
	
	while (getline(FASTQ1, line , '\n' )) 
	{	
		if(lineCt==1 )
		{
			for(int i=0; i<longestLine-merLength; ++i)
			{
				emptyVector.push_back(0);
			}
	
			for(int i=0; i<longestLine; ++i)
			{
				emptyVector2.push_back(0);
			}
			
			baseMap["A"]=emptyVector2;
			baseMap["C"]=emptyVector2;
			baseMap["G"]=emptyVector2;
			baseMap["T"]=emptyVector2;
			baseMap["N"]=emptyVector2;
			
			totalMersAtPosition=emptyVector;			
		}

		double r=rand() % 100; //Random number between 0 and 99
		
		if(lineCt%4==1 && r <samplingPercent)
		{
			map<string, int > thisReadsMers;
			
			int endPt;
			if(firstMerOnly==true)
			{
				endPt=1;
			}
			else
			{
				endPt=(int)line.size()-merLength;
			}			
			
			for(int i=0; i<endPt; ++i)
			{	
				string mer=line.substr(i, merLength);
				
				if(merMap.find(mer)==merMap.end())
				{
					merMap.insert(map< string, vector <int> >::value_type(mer, emptyVector)  );			
					merReadMap.insert(map< string, int>::value_type(mer, 0)  );		
				}
				++merMap[mer][i];	
				++totalMersAtPosition[i];
				
				//Check if it has already been inserted for this read
				if(thisReadsMers.find(mer)==thisReadsMers.end())
				{		
					thisReadsMers.insert(map< string, int>::value_type(mer, 0)  );	
					++merReadMap[mer];					
				} 
			}
			
			for(int i=0; i<line.size(); ++i)
			{
				string mer=line.substr(i, 1);
				++baseMap[mer][i];				
			}
			++totalReads;
		}
		++lineCt;
	}
	

}


void writeResults()
{
	cout<<"Writing results"<<endl;
	
	vector<int>cts;
	vector<double>totalOverRep;
	
	map< string, int > totalBase;
	map< string, double > percBase;
	
	//Determine if overRepresented 
	//lineCt is the total number if lines processed (reads*4)
	//Total expected per bin assumes random distribution
	int ctr=0;
		
	typedef std::map< string, vector <int> >::iterator it_type;
	
	//Determine if each line is over represented	
	
	vector <bool> overRep; //True false if any position in the line is over represented

	int totalBases=0;
	
	for(it_type iterator = baseMap.begin(); iterator != baseMap.end(); iterator++) 
	{
		cts=iterator->second;
		int t=std::accumulate(cts.begin(),cts.end(),0);
		totalBase[iterator->first]=t;
		totalBases=totalBases+t;
	}
	
	cout<<"total bases"<<endl;
	
	typedef std::map< string, int >::iterator it_type_int;
	
	for(it_type_int iterator = totalBase.begin(); iterator != totalBase.end(); iterator++) 
	{
		int b;
		b=iterator->second;
		percBase[iterator->first]=(double)b /(double)totalBases;
	}
	
	int nMers=emptyVector.size();	
	
	cout<<"evaluating good bad"<<endl;
	
	ctr=0;

	for(it_type iterator = merMap.begin(); iterator != merMap.end(); iterator++) 
	{
		bool bad=false;	
		double totalExpected=0;
		double totalObserved=0;
		
		vector <double> obsCt;
		
		cts=iterator->second;
		int total =std::accumulate(cts.begin(),cts.end(),0);		
		
		double expected=0;
		double percTotal=1;
	
		string thisMer=iterator->first;
		
		for(int i=0; i<thisMer.length(); ++i)
		{			
			percTotal=percTotal*percBase[thisMer.substr(i,1)];		
		}
		
		
		for(int i=0; i<cts.size(); ++i)
		{			
			double observed=(double) cts[i];	
			
			expected=percTotal*(double)totalMersAtPosition[i]; //need to calculate how many mers start at each row			
			
			obsCt.push_back(observed);
			if (((double)observed/(double)totalMersAtPosition[i]>cutOff || expected/(double)totalMersAtPosition[i]>cutOff) )			
			{
				bad=true;			
			}
			totalObserved=totalObserved+observed;
			totalExpected=totalExpected+expected;				
		}			
			
		overRep.push_back(bad);
		obsMap[thisMer]=obsCt;
		totalOverRep.push_back( totalObserved/totalExpected);	
	}
	
	cout<<"Writing results..."<<endl;
	
	//======================================================
	//Write out results
	//======================================================	
  
  
	cout<<"Writing mers..."<<endl;
	
	remove(outFileName.c_str());
	//Open output stream to write the findings for each transcript
	ofstream outputStream;
	outputStream.open(outFileName.c_str());	
	
cout<<"totalReads:"<<totalReads<<endl;
	
	ctr=0;	
	outputStream<<fastqFileName<<"\n";	
	outputStream<<"Mer\tTotalOverRep\tPercentageOfTotalReads\tPositions"<<"\n";	
	for(it_type iterator = merMap.begin(); iterator != merMap.end(); iterator++) 
	{		
		cts=iterator->second;	
		
		if(overRep[ctr]==true)
		{
			string mer=iterator->first;
			outputStream<<mer<<"\t";
			outputStream<< totalOverRep[ctr]<<"\t";
			outputStream<< (double)merReadMap[mer]/totalReads;	
				
			for(int i=0; i<cts.size(); ++i)
			{
				outputStream<< "\t"<<(double)cts[i]/(double)totalMersAtPosition[i];
			}
			outputStream<<"\n";	
		}

		++ctr;
	}
	
	outputStream.close();
	
	string outFileNameSingleBases=outFileName+".single.txt";
	remove(outFileNameSingleBases.c_str());
	//Open output stream to write the findings for each transcript
	ofstream outputStreamBase;

	outputStreamBase.open(outFileNameSingleBases.c_str());	
	
	outputStreamBase<<"Base"<< "\t"<<"Positions"<<"\n";	
	
	for(it_type iterator = baseMap.begin(); iterator != baseMap.end(); iterator++) 
	{
		outputStreamBase<< iterator->first;
		
		cts=iterator->second;
			
		for(int i=0; i<cts.size(); ++i)
		{
			outputStreamBase<< "\t"<<cts[i];
		}
		outputStreamBase<<"\n";		
		++ctr;
	}
	
	outputStreamBase.close();
		
	
}

/*===========================================================================================
Check that all of the necessary fields exist.
===========================================================================================*/

unsigned int checkErrors()
{
	//Errors 

	int err=0;
	Handy h(0);
	string problems="";
	
	if(fastqFileName.length()==0)
	{
		problems.append("A fastq file containing the reads is needed (-fastq).\n");
		++err;
	}
	if(outFileName.length()==0)
	{
		problems.append("A output file name for you bam is needed (-out).\n");
		++err;
	}
			
	/*============================================================
	Check that all the relevant files can be read/written to
	==============================================================*/

	return err;
	
}

void displayHelp()
{

	cout<<"Required parameters:\n";	
	cout<<"-fastq Name of the fastq file\n";
	cout<<"-out Name of the output file"<<endl;
	cout<<"-mer The length the KMer"<<endl;
	cout<<"Optional parameters:\n";	
	cout<<"-cutOff Only print KMers with at least (cutoff) % of reads in one cell. Default = 0.0001 "<<endl;

}

void somethingsGoneWrong(string whatsGoneWrong)
{

	cout<<"ERROR: Something has gone horribly wrong.\n";	
	cout<<whatsGoneWrong;
	cout<<"\n";	
	cerr<<"\nPlease try again.";
	
}

/*
			
			if(observed>cutOff || expected >cutOff)
			{					
				double notObserved=(double)totalMersAtPosition[i]-observed;
				double notExpected=(double)totalMersAtPosition[i]-expected;				
				
				double chi=((observed-expected)*(observed-expected) )/expected + ((notObserved-notExpected)*(notObserved-notExpected) )/notExpected;
				
				if( chi > 6.635	) //critical value of chi square with 1 df
				{
					bad=true;
				}
				
				if(chi>maxOverRep)
				{
					maxOverRep=chi;
				}
			}
			*/
