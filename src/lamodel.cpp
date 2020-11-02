// Version: $Id: lamodel.cpp 172 2014-02-12 10:06:07Z gk $
/* 
 
 
lamodel is a network simulator that simulates memory engram formation 
in a population consisting of excitatory and inhibitory neurons with independent
dendritic subunits. 

This is the entry point for the network simulator. 
This file parses the command line parameters 

Basic command line options for the simulator:

Example: ./lamodel -P 2 -T 1400 [-G] [-L] [-s datadir-name] [-S random-seed]

-P <patterns>: number of memories  to enoode
-T <minutes>: Interval between memories (minutes)
-G: Global-only protein synthesis 
-L: Dendritic-only protein synthesis 
-s <dirname>: Name of the directory to store the data (inside the data/ folder)
-S <random-seed>: Intitialize the random generator

-o <param-name>=<param-value>: Set various simulation parameters:

-o nlTypePV=[0,1,2,3]   : Set basket cell dendrite nonlinearity type. 0=supra, 1=sub, 2=linear, 3=mixed supra/sub
-o nlTypeSOM=[0,1,2,3]  : Set SOM+ cell dendrite nonlinearity type. 0=supra, 1=sub, 2=linear, 3=mixed supra/sub
-o INClustered=[0,1]    : Set whether IN synapses  should target all dendrites randomly (dispersed) or only 33% of them (clustered)


*/


#include "constructs.h"
#include <iostream>
#include <sstream>
#include <cstring>
#include <string>
#include <unistd.h>
#include <stdio.h>
#include <getopt.h>

inline float MycaDP(float x)
{
	//return x >0.2 ? 1.0 : 0.0;
	//float f =  (2.0/(1.+exp(-(x*10.-3.5)*10.))) - (1.0/(1.+exp(-(x*10.0-0.5)*19.)));


	float f = (1.3/(1.+exp(-(x*10.-3.5)*10.))) - (0.3/(1.+exp(-(x*10.0-2.0)*19.)));
	return f;
	//return f;//(2.0/(1.+exp(-(x*10.-0.7)*10.))) - (1.0/(1.+exp(-(x*10.0-0.2)*19.)));
	//return (2.0/(1.+exp(-(x*10.-3.1)*8.))) - (1.0/(1.+exp(-(x*10.0-2.1)*6.)));
}




inline void program_input2(nrn_list lst, int tstart, int duration, float freq, float randomness, int limitActive = -1)
{
	int skip = 0;
	if (limitActive == -999)
		skip = 3;
	for (nrn_iter n = lst.begin(); n != lst.end(); ++n)
	{
		if (skip>0)
		{
			skip--;
			continue;
		}

		LAInput* in = (LAInput*)(*n);
		in->Program(tstart, duration, freq, randomness);
		//if (limitActive >0 && ++tot >= limitActive) return;
	}
}



void RunPattern2(LANetwork& net,  int pat[])
{

	int pad = 100;
	int duration = 3800;

	int i;

	for (i=0; i < 100; i++)
	{
		LAInput* in = (LAInput*) net.pattern_inputs[i];
		if (pat[i]>0)
			in->Program(pad, duration, 80, 0.5);
		else
			in->Program(pad, duration, 0, 0.5);
	}

	program_input2(net.bg_list, pad, duration, .5, 0.5, -1);


	net.ResetSpikeCounters();

	for (syn_iter si = net.synapses.begin(); si!= net.synapses.end(); si++)
		(*si)->calcium =0.;

	net.StimDynamics(duration+pad+pad);
	
	int tActive =0;
	int tSpikes =0;
	
	for (nrn_iter ni = net.pyr_list.begin(); ni != net.pyr_list.end(); ni++)
	{
		LANeuron* n = *ni;
		if (float(n->total_spikes) /4. > 5. )
			tActive++;
		tSpikes += n->total_spikes;
	}

	net.totPyrActive = tActive;
	
	printf("Active pyrs= %d (%.2f%%), mean ff= %.2f\n", tActive, 100.0*float(tActive)/float(net.pyr_list.size()), tSpikes/(float(net.pyr_list.size())*2));
}


void MyInterstim(LANetwork& net, int durationSecs)
{

	float totWeightUp = 0.0, totWeightDown =0.0;

	// Generate synaptic tags
	//
	cout << "Neuron calc: ";
	for (nrn_iter ni = net.pyr_list.begin(); ni != net.pyr_list.end(); ni++)
	{
		LANeuron*  nrn = *ni;
		float nrnCalc =0.0;

		if (nrn->type  != 'P') continue;


		for (branch_iter bi = nrn->branches.begin(); bi != nrn->branches.end(); ++bi)
		{
			LABranch* b = *bi;
			for (syn_iter si = b->synapses.begin(); si != b->synapses.end(); ++si)
			{
				LASynapse* s =*si;
				b->totcalc += s->calcium;
			}

			nrnCalc +=  b->totcalc;

		}



		nrn->totcalc  = nrnCalc;

		//cout  << nrn->total_spikes << " ("<< nrnCalc << ") ";

		if (nrn->totcalc > 18.0)
		{

			for (branch_iter b = nrn->branches.begin(); b != nrn->branches.end(); ++b)
				for (syn_iter si = (*b)->synapses.begin(); si != (*b)->synapses.end(); ++si)
				{
					LASynapse* s =*si;
					float ctag = MycaDP(s->calcium/40.);

					if (ctag > 0.1 || ctag < -0.1)
					{
						if (rgen()   < 0.5)
						{
							if (ctag>0) 
							{
								totWeightUp += (0.8 - s->weight);
								s->weight = 0.8;
							}
							else 
							{
								totWeightDown += (s->weight);
								s->weight = 0.1;
							}
						}
					}

				}
		}


	}
	cout <<endl;


	printf("Tot went up: %f, down: %f\n", totWeightUp, totWeightDown);

}


template <typename T> static void PrintVector2( vector<T>&  ar, ostream& outfile) 
{
	for (typename vector<T>::iterator it = ar.begin(); it != ar.end(); it++)
	{
		outfile << *it << ' ';
	}
	outfile << std::endl;
}





void RunIteration(LANetwork& net, int run_no)
{

	//int bpat[TOT_PAT][100];
	unsigned int i;

	string dirname = "./";

	ifstream inpFile("inspikes.txt");

	vector < vector <int> > spikeTimings;

	int row = 0;
	int  val;
	while(!inpFile.eof())
	{
	    std::string str = "";
	    std::getline(inpFile, str);

	    if(str[0] == ';') continue; //for comments

		if (str.length() < 2 ) continue;

	    std::stringstream ss(str);
	    int tspike=0;

		vector<int> spikes;
	    while(ss >> val) 
		{
			if (val >0) 
				spikes.push_back(tspike);

			tspike++;
		}
	    spikeTimings.push_back( spikes);
		row++;
	}

	cout << "Input spikes: " <<endl;
	for (i=0; i < spikeTimings.size(); i++)
	{
		//for (k=0; k < spikeTimings.at(i).size(); k++)
		//	cout << spikeTimings.at(i).at(k) << " " ;
		cout << spikeTimings.at(i).size() << " ";
	}
	cout<<endl;


	net.enablePlasticity = true;
	net.ResetSpikeCounters();

	for (syn_iter si = net.synapses.begin(); si!= net.synapses.end(); si++)
		(*si)->calcium =0.;

	
	for (i=0; i < spikeTimings.size() && i <100; i++)
	{
		LAInput* in = (LAInput*) net.pattern_inputs.at(i);
		in->SetSpikes(spikeTimings.at(i));
	}

	cout << "Running training (10000msec)..." << endl;
	net.StimDynamics(10000);

	cout <<"Output spikes: " ;
	for (nrn_iter ni = net.pyr_list.begin(); ni != net.pyr_list.end(); ni++)
	{
		LANeuron* n = *ni;
		cout << n->total_spikes << " ";
	}
	cout <<endl;

	cout << "Saving spikes.txt..." ;
	ofstream spikesdat((dirname + "spikes.txt").c_str());

	for (nrn_iter na = net.pyr_list.begin(); na != net.pyr_list.end(); ++na)
	{
		LANeuron* nrn = *na;

		int spikes[10000] = {0};
		for (vector<int>::iterator ii = nrn->spikeTimings.begin(); ii != nrn->spikeTimings.end(); ii++)
		{
			int t = *ii;
			if (t>0 && t < 10000)
				spikes[t] = 1;
		}

		for (i=0; i < 10000; i++)
			spikesdat << spikes[i] << " ";
		spikesdat << endl;
	}

	spikesdat.close();
	cout << "Done." <<endl;


	cout << "Interstim (3600)"  << endl ;
	MyInterstim(net, (3600));
	cout << "Done."<<endl;

/*
	ofstream traincnt(("train_spikes_"+suffix+".txt").c_str());
	ofstream caapcnt(("caap_spikes_"+suffix+".txt").c_str());
	ofstream pyrActiveCnt(("pyr_active"+suffix+".txt").c_str());

	for (i=0; i < TOT_PAT; i++)
	{
		net.ResetSpikeCounters();

		cout << "Running " << i << " "   ;
		for (k=0; k < 100; k++) cout <<bpat[i][k] ; 
		cout <<endl;

		RunPattern2(net, bpat[i]);
		for (nrn_iter i = net.neurons.begin(); i != net.neurons.end(); i++)
		{
			traincnt<< (*i)->total_spikes << " ";
			//caapcnt<< (*i)->totalCaap << " ";
		}
		traincnt <<endl;

		pyrActiveCnt << net.totPyrActive << endl;



		cout << "Interstim "  << endl ;
		MyInterstim(net, (1000));
	}

	net.enablePlasticity = false;


	ofstream testcnt(("test_spikes"+suffix+".txt").c_str());
	for (i=0; i < TOT_PAT; i++)
	{
		net.ResetSpikeCounters();
		cout << i << endl;
		RunPattern2(net,bpat[i]);
		
		for (nrn_iter i = net.neurons.begin(); i != net.neurons.end(); i++)
		{
			testcnt<< (*i)->total_spikes << " " ;
			//cout << (*i)->total_spikes << " " ;
		}
		testcnt <<endl;
		cout << endl ;
		

		pyrActiveCnt << net.totPyrActive << endl;
	}
	*/
}





// Parse command line  arguments
int main( int argc, char* argv[])
{
	
	int c;

	// Set defaults
	int nneurons = 500;
	int nbranches = 20;
	int ninputs = 400;
	int nperinput = 1;
	int npatterns = ninputs;
	int nonesperpattern = 1;
	int interstim = 60;
	int rseed = 1980;
	char* suffix = NULL;
	bool disableCreb = false;
	int numPatterns = 4;

	LANetwork net; 
	net.enablePruning = true; // default

	while ((c = getopt(argc, argv, "M:N:H:B:I:i:P:p:T:S:s:d:w:O:g:l:b:c:o:t:xnLDRJCGhU"))!= -1)
	{
		switch (c)
		{
			case '?':
			case 'h':
			cout <<
			"Basic command line options for the simulator:" << endl << 
			"Example: ./lamodel -P 2 -T 1400 [-G] [-L] [-s datadir-name] [-S random-seed]" << endl << 
			"-P <patterns>: number of memories  to enoode"<< endl << 
			"-T <minutes>: Interval between memories (minutes)"<< endl << 
			"-G: Global-only protein synthesis "<< endl << 
			"-L: Dendritic-only protein synthesis "<< endl << 
			"-s <dirname>: Name of the directory to store the data (inside the data/ folder)"<< endl << 
			"-S <random-seed>: Intitialize the random generator"<< endl << 
			""<< endl << 
			"-o <param-name>=<param-value>: Set various simulation parameters:"<< endl << 
			""<< endl << 
			"-o nlTypePV=[0,1,2,3]   : Set basket cell dendrite nonlinearity type. 0=supra, 1=sub, 2=linear, 3=mixed supra/sub"<< endl << 
			"-o nlTypeSOM=[0,1,2,3]  : Set SOM+ cell dendrite nonlinearity type. 0=supra, 1=sub, 2=linear, 3=mixed supra/sub" << endl << 
			"-o INClustered=[0,1]    : Set whether IN synapses  should target all dendrites randomly (dispersed) or only 33% of them (clustered)" << endl;
			return 1;
			break;

			case 'B': nbranches = atoi(optarg); break;
			//case 'I': ninputs = atoi(optarg); break;
			//case 'i': nperinput = atoi(optarg); break;
			case 'N': nneurons = atoi(optarg); break;
			case 'P': npatterns = atoi(optarg); break;
			case 'p': nonesperpattern = atoi(optarg); break;
			case 'T': interstim = atoi(optarg); break;
			case 'S': rseed = ( atoi(optarg)); break;
			case 's': suffix = strdup(optarg); break;

			case 'n': disableCreb = true; break;
			case 'w': net.isWeakMem.push_back(atoi(optarg)-1); break;

			case 'L': net.localProteins = true; break;
			case 'G': net.globalProteins = true; break;
			case 'D': net.debugMode = true; break;
			case 'R': net.repeatedLearning = true; break;
			case 'J': net.pretraining = true; break;
			case 'C': net.altConnectivity = true; break;
			case 'O': net.branchOverlap = atof(optarg); break;
			case 'H': net.homeostasisTime = atof(optarg); break;


			case 'o': 
				char* o = strstr(optarg, "=");
				if (o)
				{
					*o = '\0';
					char* val = o+1;

					if (!strcmp(optarg, "connectivityParam")) net.connectivityParam = atof(val); 
					else if (!strcmp(optarg,  "BSPTimeParam")) net.BSPTimeParam = atof(val); 
					else if (!strcmp(optarg,  "homeostasisTimeParam")) net.homeostasisTimeParam = atof(val); 
					else if (!strcmp(optarg,  "CREBTimeParam")) net.CREBTimeParam = atof(val); 
					else if (!strcmp(optarg,  "inhibitionParam")) net.inhibitionParam = atof(val); 
					else if (!strcmp(optarg,  "globalPRPThresh")) net.globalPRPThresh = atof(val); 
					else if (!strcmp(optarg,  "localPRPThresh")) net.localPRPThresh = atof(val); 
					else if (!strcmp(optarg,  "dendSpikeThresh")) net.dendSpikeThresh = atof(val); 
					else if (!strcmp(optarg,  "initWeight")) net.initWeight*= atof(val); 
					else if (!strcmp(optarg,  "maxWeight")) net.maxWeight*= atof(val); 
					else if (!strcmp(optarg,  "stimDurationParam")) net.stimDurationParam = atof(val); 
					else if (!strcmp(optarg,  "nNeuronsParam")) nneurons *= atof(val); 
					else if (!strcmp(optarg,  "nBranchesParam")) nbranches *= atof(val); 
					else if (!strcmp(optarg,  "nBranchesTurnover")) net.nBranchesTurnover = atoi(val); 
					else if (!strcmp(optarg,  "INClustered")) net.INClustered = atoi(val); 
					else if (!strcmp(optarg,  "nlTypePV")) net.nlTypePV = atoi(val); 
					else if (!strcmp(optarg,  "nlTypeSOM")) net.nlTypeSOM = atoi(val); 
					else if (!strcmp(optarg,  "hasCaap")) net.hasCaap = atoi(val); 
					else if (!strcmp(optarg,  "numPatterns")) numPatterns = atoi(val); 

					printf("Param name='%s' value='%f'\n", optarg, atof(val));
				}
			break;
		}
	}

	LANetwork::SetRandomSeed(rseed);
	net.disableCreb = disableCreb;

	// Create the network 
	net.CreateFearNet(nneurons, nbranches, ninputs, nperinput);

	char buf[512];
	if (suffix)
		sprintf(buf, "./data/%s", suffix );
	else
		sprintf(buf, "./data/N%d.B%d.I%d.i%d.P%d.p%d.T%d.S%d.w%d_%s", nneurons, nbranches, ninputs, nperinput, npatterns, nonesperpattern, interstim, rseed, (int)net.isWeakMem.size(),  suffix ? suffix : "");

	cout << "Data directory is "<< buf <<  endl;
	net.SetDataDir( buf);


	struct stat  mstat;
	int itNo =0;
	time_t lastT =0;

	if (stat("inspikes.txt", &mstat) != -1) 
		lastT = mstat.st_mtime;

	cout << "Waiting for file 'inspikes.txt' to be updated, press Ctrl-C to exit"<<endl;
	while (true)
	{
		if (stat("inspikes.txt", &mstat) == -1) 
		{
			cout << "inspikes.txt not found" << endl; 
		}
		else if (lastT != mstat.st_mtime)
		{
			lastT = mstat.st_mtime;
			printf("[%d] Running iteration...\n", itNo);

			RunIteration(net, itNo);
			printf("[%d] Done. Waiting for file...\n", itNo);
			itNo++;
		}
		sleep(1);
	}

	return 0;

}



