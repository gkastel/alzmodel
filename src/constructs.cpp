/* 

lamodel main implementation file. 
See constructs.h for class definitions

*/



#include "constructs.h"
#include <iterator>
#include <assert.h>

const int DEBUG_SID = 3024;
const int CUTOFF = 10.0;


#define VEC_REMOVE(vec, t) (vec).erase(std::remove((vec).begin(), (vec).end(), t), (vec).end())


// PROTEIN SYNTHESIS thresholds. When calcium crosses this value, proteins will be synthesized
const float GPROD_CUTOFF = 18.0; // 18.0 Global ca threshold
const float BPROD_CUTOFF = 1.8;  // 1.8 Dendrite ca threshold


int LANetwork::RSEED = 1980;



// The alpha function for protein sythesis over time (x is time)
inline double nuclearproteinalpha(float x)
{
	return (x>20.)*((x-20.*60.)/(30.*60)) * exp(1. - (x-20.*60. )/(30.*60.));
}




// The alpha function for protein sythesis in dend branches over time (x is time)
inline double branchproteinalpha(float x)
{
	return ((x)/(15.*60)) * exp(1. - (x )/(15.*60.));

	
}



// the curve for the magnitude of LTP vs Ca++  . (x is calcium)
inline float caDP(float x)
{
	//return x >0.2 ? 1.0 : 0.0;
	//float f =  (2.0/(1.+exp(-(x*10.-3.5)*10.))) - (1.0/(1.+exp(-(x*10.0-0.5)*19.)));


	float f = (1.3/(1.+exp(-(x*10.-3.5)*10.))) - (0.3/(1.+exp(-(x*10.0-2.0)*19.)));
	return f;
	//return f;//(2.0/(1.+exp(-(x*10.-0.7)*10.))) - (1.0/(1.+exp(-(x*10.0-0.2)*19.)));
	//return (2.0/(1.+exp(-(x*10.-3.1)*8.))) - (1.0/(1.+exp(-(x*10.0-2.1)*6.)));
}


/*
static void clampval(float&  cmin, float& cmax, float val)
{
	if (val < cmin )
		cmin = val;
	if (val > cmax )
		cmax = val;
}
*/






// Preallocates spikes in a list  of neurons
inline void program_input(nrn_list lst, int tstart, int duration, float freq, float randomness, int limitActive = -1)
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




// Create a list of neurons
void LANetwork::CreateNeurons(int number, int n_branches_per_neuron, char type, vector<LANeuron*>* appendTo = 0, int inputId =-1, int somethingDummy = 0)
{
	for (int i =0 ;i < number; i++)
	{
		LANeuron* n ;
		if (type == 'S')
		{
			n = new LAInput();
			n->network = this;

			n->input_id = inputId;
			((LAInput*)n)->groupIdx = i;
		}
		else
		{
			n = new LANeuron();
			n->network = this;
			for (int tt =0; tt < n_branches_per_neuron; tt++)
			{
				LABranch* bb  = new LABranch;
				bb->bid = this->branches.size();
				bb->neuron = n;


				// Set type of nonlinearity according to dendrite type and command line options

				if (type == 'P')  // Pyramidals
				{
					bb->nlType = DEND_SUPRA;
				}
				else if ( type == 'M') // SOM interneurons
				{
					if (nlTypeSOM == DEND_MIXED)
					{
						if (tt < 0.5*n_branches_per_neuron) bb->nlType = DEND_SUB;
						else bb->nlType = DEND_SUPRA;
					}
					else if (nlTypeSOM < DEND_MIXED && nlTypePV >= 0)
						bb->nlType = nlTypeSOM;
					else
					{
						printf("bad type %d", nlTypeSOM);
						abort();
					}
				}
				else if (type == 'V') // basket interneurons
				{
					if (nlTypePV == DEND_MIXED)
					{
						if (tt < 0.5*n_branches_per_neuron) bb->nlType = DEND_SUB;
						else bb->nlType = DEND_SUPRA;
					}
					else if (nlTypePV < DEND_MIXED && nlTypePV >= 0)
						bb->nlType = nlTypePV;
					else
					{
						printf("bad type %d", nlTypePV);
						abort();
					}

				}


				this->branches.push_back(bb);
				n->branches.push_back(bb);
			}
		}
		n->type = type;

		n->V = param_V_rest +1.0;
		n->nid = this->neurons.size();

		n->glx = 1.9*(n->nid%50)/50. - 0.9;
		n->gly = 1.9*(n->nid/50)/50. - 0.9;

		this->neurons.push_back(n);
		if (appendTo)
			appendTo->push_back(n);
	}
}




void LANetwork::CalculateDistances()
{
	for (nrn_iter na = neurons.begin(); na != neurons.end(); ++na)
		for (nrn_iter nb = neurons.begin(); nb != neurons.end(); ++nb)
			if (nb != na)
			{
				LANeuron* a = *na, *b=*nb;

				
				float dx = b->pos_x - a->pos_x;
				float dy = b->pos_y - a->pos_y;
				float dz = b->pos_z - a->pos_z;
				distances[pair<int,int>(a->nid,b->nid)] = (double)sqrt(dx*dx + dy*dy + dz*dz);
			}
}



void LANetwork::AddSynapse(LANeuron* a, LABranch* br, float weight, bool isPlastic)
{
		LASynapse* syn = new LASynapse();
		syn->sid  = this->synapsesCounter++;
		syn->source_nrn = a;
		syn->target_nrn = br->neuron;
		syn->isPlastic = isPlastic;

		syn->weight  = weight;
		syn->target_branch = br; 

		syn->source_nrn->outgoing.push_back(syn);
		syn->target_nrn->incoming.push_back(syn);
		syn->target_branch->synapses.push_back(syn);
		syn->pos = rgen();
		this->synapses.push_back(syn);
}



// Connect two lists of neurons
// Distances not implemented in this version
int LANetwork::ConnectNeurons(vector<LANeuron*> fromList, vector<LANeuron*> toList, bool isClustered, float toDistance, int nNeuronPairs, int nSynapsesPerNeuron, float weight, bool isPlastic= false, bool randomizeweight = false, float overlap =-1.0)
{
	int tpairs =0;
	while(true)
	{
		LANeuron* a = fromList.at(int(rgen()*(float)fromList.size()));
		LANeuron* b = toList.at(int(rgen()*(float)toList.size()));

		for (int i =0; i < nSynapsesPerNeuron; i++)
		{

			float rval;
			if (isClustered)
				rval = rgen()*float(b->branches.size()/3);  // Clusterd
			else
				rval = rgen()*float(b->branches.size());

			LABranch* br =  b->branches[(int)rval];

			this->AddSynapse(a, br, weight, isPlastic);	
		}

		if (++tpairs >= nNeuronPairs) break;
	}
	return tpairs;
}




inline int randomindex(int max)
{
	return (int)(rgen()*float(max));
}



// Remove input synapses (turnover)
int LANetwork::PurgeInputSynapses(int totalToRemove,float weightLimit)
{

	int totalFound =0;
	int totalTries = totalToRemove*3;
	while (totalFound < totalToRemove && totalTries-->0)
	{
		nrn_list lst = this->inputs_cs.at(randomindex(this->inputs_cs.size()));
		LANeuron* n = lst.at(randomindex(lst.size()));
		if (n->outgoing.size())
		{
			LASynapse* s = n->outgoing.at(randomindex(n->outgoing.size()));
			if (s->target_nrn->type == 'P'  && s->weight <= weightLimit)
			{
				//candidate for deletion

				//pthread_mutex_lock(&this->synapses_mutex);

				VEC_REMOVE(s->source_nrn->outgoing, s);
				VEC_REMOVE(s->target_nrn->incoming, s);
				VEC_REMOVE(s->target_branch->synapses, s);
				VEC_REMOVE(this->synapses, s);
				//cout << " Removing " << s->sid << endl;
				delete s; 

				//pthread_mutex_unlock(&this->synapses_mutex);

				totalFound++;
			}
		}
	}

	if (totalFound < totalToRemove)
	{
		cout << " Warning: not enougn synapses to remove: " << totalToRemove << endl;  
	}

	return totalFound;
}





// Construct the network
void LANetwork::CreateFearNet(int nneurons, int nbranches, int ninputs, int nneuronsperinput)
{
	this->n_neurons = nneurons;
	this->n_inputs = ninputs;
	this->n_neurons_per_input = nneuronsperinput;
	this->n_branches_per_neuron = nbranches;

	// Excitatory population (Pyr)
	this->CreateNeurons(this->n_neurons*0.8, this->n_branches_per_neuron, 'P', &this->pyr_list, -1, this->nBranchesTurnover);

	// Inhibitory populations 
	this->CreateNeurons(this->n_neurons*0.1, 10 , 'V', &this->in_pv); // PV

	this->CreateNeurons(this->n_neurons*0.1, 10 , 'M', &this->in_som); // SOM


	
	CreateNeurons(100, 0, 'S', &this->pattern_inputs, -1, true);


	// Background noise inputs (not used)
	CreateNeurons(10, 0, 'S', &this->bg_list, -1, true);

	this->CalculateDistances();

	float  baseSyns = 1000;
	baseSyns = 500;


	// Pur <-> basket cells
	ConnectNeurons(this->pyr_list, this->in_pv, this->INClustered, 10.,  1*baseSyns, 1, .3, true);
	ConnectNeurons(this->in_pv, this->pyr_list, 0, 10.,  1*10.*baseSyns, 1, .3, false);

	// Pyr <-> SOM cells
	ConnectNeurons(this->pyr_list, this->in_som, 0, 10.,  2*baseSyns, 1, .3, true);
	ConnectNeurons(this->in_som, this->pyr_list, 0, 10., 4*baseSyns, 1, .3, true);

	//baseSyns = 8000;
	baseSyns = 3500;


	// Memory inputs -> Pyr 
	this->ConnectNeurons(this->pattern_inputs, this->pyr_list, 0, 10.0, baseSyns, 1, 0.3 /* initWeight */, true, false, false);


	this->RecordInitialWeightSums();
}










template <typename T> static void PrintVector( vector<T>&  ar, ostream& outfile) 
{
	for (typename vector<T>::iterator it = ar.begin(); it != ar.end(); it++)
	{
		outfile << *it << ' ';
	}
	outfile << std::endl;
}




void LANetwork::StoreDataFiles( bool extras = false )
{
	vector<int> pattern;

	string dirname = this->datadir;
	ofstream paramsdat((dirname + "/parameters.txt").c_str());
	paramsdat <<"total_neurons="<< this->neurons.size() << endl;
	paramsdat <<"total_pyramidals="<< this->pyr_list.size() << endl;
	paramsdat <<"branches_per_neuron="<< this->n_branches_per_neuron << endl;
	paramsdat <<"number_inputs="<< this->inputs_cs.size() << endl;
	paramsdat <<"neurons_per_input="<< this->n_neurons_per_input << endl;
	paramsdat <<"rseed="<< RSEED << endl;


	ofstream patternsdat((dirname + "/patterns.txt").c_str());
	for (vector<vector<int> >::iterator it = this->patterns.begin(); it != this->patterns.end(); it++)
	{
		pattern = *it;
		copy(pattern.begin(), pattern.end(), ostream_iterator<int>(patternsdat, " "));
		patternsdat << endl;
	}

	ofstream synstatedat((dirname + "/synstate.dat").c_str());

	ofstream spikesdat((dirname + "/spikes.dat").c_str());
	ofstream crebdat((dirname + "/creb.dat").c_str());
	ofstream voltagedat((dirname + "/voltages.dat").c_str());
	ofstream branchspikesdat((dirname +"/branchspikes.dat").c_str());
	ofstream branchcalcium((dirname + "/branchcalcium.dat").c_str());
	ofstream weightsdat((dirname + "/weights.dat").c_str());
	ofstream branchproteins((dirname + "/branchproteins.dat").c_str());
	ofstream branchstrengths((dirname + "/branchstrengths.dat").c_str());
	ofstream tagsdat((dirname + "/tags.dat").c_str());
	ofstream nrnproteindat((dirname + "/nrnprotein.dat").c_str());
	ofstream weighthistorydat((dirname + "/weighthistory.dat").c_str());
	ofstream dbgneuron((dirname + "/dbgneuron.dat").c_str());

	ofstream syn_per_branch((dirname + "/syn_per_branch.dat").c_str());

	PrintVector<float>( dbgNeuron, dbgneuron);

	for (nrn_iter na = neurons.begin(); na != neurons.end(); ++na)
	{
		LANeuron* nrn = *na;

		PrintVector<int>( nrn->spikeTimings, spikesdat);
		PrintVector<float>( nrn->proteinHistory, nrnproteindat);

		for (branch_iter bi = nrn->branches.begin(); bi != nrn->branches.end(); ++bi)
		{
			LABranch* b = *bi;

			PrintVector<float>(b->branchSpikesHistory, branchspikesdat);

			for (syn_iter si = b->synapses.begin(); si != b->synapses.end(); ++si)
			{
				LASynapse* s = *si;
				//PrintVector<float>(s->weightHistory, weightsdat);
				synstatedat << s->sid<<" " 
					<< b->bid<<" "
					<< nrn->nid  << " "
					<< s->source_nrn->nid <<" " 
					<< s->source_nrn->input_id<< " "
					<< b->strength  << " "
					<< s->weight << " " 
					<<endl;

				if (s->tagHistory.size())
				{
					tagsdat << s->source_nrn->input_id<< " ";
					//PrintVector<float>(s->tagHistory, tagsdat);
				}

				if (s->weightHistory.size())
				{
					weighthistorydat << s->source_nrn->input_id<< " ";
					//PrintVector<float>(s->weightHistory, weighthistorydat);
				}

				if (s->isPlastic && s->source_nrn->input_id >=0 && s->weight > .7)
				{
				//	totPot[ s->source_nrn->input_id >=0 ] += 1;
				}
			}

		}
	}

	ofstream sppdat((dirname + "/spikesperpattern.dat").c_str());
	for (uint i =0; i < this->spikesPerPattern.size(); i++)
	{
		for (uint j =0; j < this->spikesPerPattern[i].size(); j++) sppdat << this->spikesPerPattern[i][j] << " ";
		sppdat << endl;
	}



}


// Stimulation dynamics, with dt=1msec 
void LANetwork::StimDynamics(int duration) 
{
	int t = 0;
	bool spikeState[this->neurons.size()+1];
	int lastSpikeT[this->neurons.size()+1];


	fill_n(spikeState, this->neurons.size()+1, 0);
	fill_n(lastSpikeT, this->neurons.size()+1, 0);


	for (syn_iter si = this->synapses.begin(); si != this->synapses.end(); ++si)
	{
		(*si)->calcium = 0.0;
	}


	for (nrn_iter ni = this->neurons.begin(); ni != this->neurons.end(); ++ni)
	{
		LANeuron* n = *ni;
		n->wadapt = 0.0;
		n->vspike =0.;
		n->vreset =0.;
		n->V =0.;
	}



	for (branch_iter bi=this->branches.begin(); bi != this->branches.end(); ++bi)
	{
		(*bi)->totcalc = 0.0;
		(*bi)->depol = 0.0;
		(*bi)->dspike = 0.0;
		(*bi)->dspikeT = -1;
	}

	for (t=0; t < duration; t++)
	{
		for (nrn_iter ni = this->neurons.begin(); ni != this->neurons.end(); ++ni)
		{
			LANeuron* n = *ni;
			float soma_inh =0;
			float soma_exc =0;

			for (branch_iter bi=n->branches.begin(); bi != n->branches.end(); ++bi)
			{
				LABranch* b = *bi;
				float dend_exc =0.;
				float dend_inh =0.;
				for (syn_iter si=b->synapses.begin(); si != b->synapses.end(); ++si)
				{
					LASynapse* s = *si;
					if (spikeState[s->source_nrn->nid])
					{
						if (s->source_nrn->type == 'V' ) dend_inh += ( s->weight);
						else if (s->source_nrn->type == 'M' ) soma_inh += ( s->weight);
						else dend_exc += (  s->weight);
					}
				}

				if (b->nlType == DEND_SUB)
				{
					// sublinear integration 
					b->depol +=  pow(4.0*dend_exc - 3.*dend_inh, 0.7) -  b->depol/20.0;
				}
				else if (b->nlType == DEND_SUPRA) 
				{
					// supralinear integration 
					b->depol +=  (4.0*dend_exc - 3. * dend_inh) -  b->depol/20.0; 
					if (b->dspikeT < t-70 && (n->vspike + b->depol) > this->dendSpikeThresh) // Generate a dendritic branch spike
					{
						b->depol = 50;
						b->dspikeT =t;
						b->branch_spikes++;
						n->dend_spikes++;
					}
				}
				else
				{
					// linear integration 
					b->depol +=  (4.0*dend_exc - 3.*dend_inh) -  b->depol/20.0;
				}

				// sum up calcium influx
				for (syn_iter si=b->synapses.begin(); si != b->synapses.end(); ++si)
				{
					LASynapse* s = *si;
					if (spikeState[s->source_nrn->nid])
					{
						if (this->enablePlasticity && s->isPlastic) 
						{
							float depol =  b->depol + n->vspike;
							if (depol > 1.0)
							{
								float ff =  (1.0/(1.+exp( (-(depol-30.0)/5.0))));
								s->calcium +=  1.1*ff;
							}
						}
					}
				}
				soma_exc +=  (b->depol); 
			}

			soma_exc *= 0.12;
			soma_inh *= 0.18;


			if (n->type == 'S') // Source neuron
			{
				LAInput* in = (LAInput*)n;
				if (in->spikeTimes.size()>0 && t >= in->nextSpikeT && in->nextSpikeT >0)
				{
					// Emit a spike
					if (in->curSpike < in->totalSpikes)
						in->nextSpikeT = in->spikeTimes.at(in->curSpike++);
					else 
						in->nextSpikeT = -1;

					spikeState[in->nid] = 1;
					n->isSpiking = true;
					lastSpikeT[in->nid] = t;
					in->total_spikes++;
					in->V = 20.0; // threshold

				}
				else
				{
					spikeState[in->nid] = 0;
					n->isSpiking = false;
					in->V = param_E_L;
				}
			}
			else /// Pyr or interneuron
			{


				if (spikeState[n->nid])
				{
					// Voltage reset
					n->V = 0.0;
					n->wadapt += 0.18;
				}

				if (n->type == 'V') // PV interneuron
	 			{
					n->V +=  (soma_exc - soma_inh) - (n->V)/10.; // - n->wadapt*(n->V+10.0);
				}
				else if (n->type == 'M') // SOM interneuron
	 			{
					n->V +=  (soma_exc - soma_inh) - (n->V)/10.; // - n->wadapt*(n->V+10.0);
				}
				else
				{
					// Pyr neuron
					//

					// dCap filter				
					//
					
					if (this->hasCaap && n->type == 'P') 
					{
						if (t - n->dCaapT < 20)
						{
							//float d = t - n->dCaapT;
							//soma_exc = 2. + 3.5 *d/20.;
							soma_exc = 3.0;
						} 
						else if (n->dCaapT < t-200 && (soma_exc > 2.5 && soma_exc < 3.) ) 
						{

							n->dCaapT = t;
							n->totalCaap ++;

						}

					}

					


					n->V +=  soma_exc - 3.0*soma_inh - (n->V)/30. -   n->wadapt*(n->V+10.0) ;

					if (this->disableCreb)
						n->wadapt -= n->wadapt/180.;
					else
						n->wadapt -= n->wadapt/((180. - 70.0*(n->crebLevel>0.2 ? 1. : 0.)));
				}


				if ( lastSpikeT[n->nid] < t-2 && n->V > (20.0 - (n->crebLevel>100. ? 2.0 : 0.) ))
				{
					// Generate a spike
					spikeState[n->nid] = 1;
					lastSpikeT[n->nid] = t;
					n->total_spikes++;
					n->isSpiking = true;
					n->vspike = 30.0; 
					n->V = 70.;

				}
				else
				{
					if (n->isSpiking)
					{
						spikeState[n->nid] = 0;
						n->isSpiking = false;
					}
				}

				// backpropagating spike
				if (n->vspike > 0) n->vspike -= n->vspike / 17.0; 

				// adaptation parameter
				if (n->wadapt <-0.) n->wadapt =0.;

			}

			if (spikeState[n->nid])
			{
				// Record this spike
				n->spikeTimings.push_back(t+ this->Tstimulation);
			}
		}

		#ifdef WXGLMODEL_H
		if (this->wx_window && t%10==0)
		{
			this->wx_window->UpdatePanel();
		}
		#endif
		

	}


	this->Tstimulation += duration;
}




void LANetwork::SaveSnapshot(char* filename)
{
	ofstream of(filename);
	for (nrn_iter ni = this->pyr_list.begin(); ni != this->pyr_list.end(); ni++)
	{
		LANeuron* nrn = *ni;
		if (nrn->type == 'P')
		{
			float totTag =0;
			float totProtein =0;
			for (branch_iter bi = nrn->branches.begin(); bi != nrn->branches.end(); ++bi)
			{
				LABranch* b = *bi;
				for (syn_iter si = b->synapses.begin(); si != b->synapses.end(); ++si)
				{
					LASynapse* s = *si;
					if (s->isPlastic)
					{
						totTag += s->stag;
					}
				}
				totProtein += b->protein;
			}
			of << nrn->nid << ' ' << totProtein << ' '<<  totTag << endl;
		}
	}

}



void LANetwork::Interstim(int durationSecs)
{

	float totWeightUp = 0.0, totWeightDown =0.0;

	// Generate synaptic tags
	//
	cout << "Neuron calc: ";
	for (nrn_iter ni = this->pyr_list.begin(); ni != this->pyr_list.end(); ni++)
	{
		LANeuron*  nrn = *ni;
		float nrnCalc =0.0;


		if (nrn->type  == 'P')
		for (branch_iter bi = nrn->branches.begin(); bi != nrn->branches.end(); ++bi)
		{
			LABranch* b = *bi;

			if (!this->enablePlasticity)
				continue;

			for (syn_iter si = b->synapses.begin(); si != b->synapses.end(); ++si)
			{
				LASynapse* s =*si;
				b->totcalc += s->calcium;
			}


			
			// Branch PRPs 
			if ( b->totcalc  > this->localPRPThresh*BPROD_CUTOFF) // This branch should produce PRPs now 
			{
				b->prpTransients.push_back( pair<float,float>(T, b->totcalc));
			}

			nrnCalc +=  b->totcalc;

		}



		nrn->totcalc  = nrnCalc;

		cout  << nrn->total_spikes << " ("<< nrnCalc << ") ";

		if (nrn->totcalc > 14.0)
		{

			for (branch_iter b = nrn->branches.begin(); b != nrn->branches.end(); ++b)
				for (syn_iter si = (*b)->synapses.begin(); si != (*b)->synapses.end(); ++si)
				{
					LASynapse* s =*si;
					float ctag = caDP(s->calcium);

					if (ctag > 0.1 || ctag < -0.1)
					{
						if (rgen()   < 0.5)
						{
							s->weight = ctag;

							cout << ctag;
							if (ctag>0.) totWeightUp += ctag;
							else  totWeightDown += ctag;
						}
					}

				}
		}


	}
	cout <<endl;


	printf("Tot went up: %f, down: %f\n", totWeightUp, totWeightDown);

}





// Perform synapse turnover
void LANetwork::DoTurnover(float durationSecs )
{

	float pTurnOver = 4.*durationSecs/86400.; // per day

	int synsAdded=0;
	for (branch_iter bi = this->branches.begin(); bi != this->branches.end(); bi++)
	{
		LABranch* b = *bi;
		if (b->turnoverRate>0.0)
		{
			vector<LASynapse*> v;
			for (syn_iter si = b->synapses.begin(); si != b->synapses.end(); si++)
			{
				LASynapse* s = *si;
				if (s->isPlastic && s->weight <= this->initWeight)
				{
					float p;
					float dist = s->pos - b->turnoverPoint;
					p = exp(-(dist*dist)/.1)*b->turnoverRate;

					if (rgen() < p*pTurnOver)
					{
						VEC_REMOVE(s->source_nrn->outgoing, s);
						VEC_REMOVE(s->target_nrn->incoming, s);
						VEC_REMOVE(this->synapses, s);

						si = b->synapses.erase(si); //VEC_REMOVE(s->target_branch->synapses, s);

						delete s; 

						/* Connect  a random input back here */
						nrn_list lst = this->inputs_cs.at(randomindex(this->inputs_cs.size()));
						LANeuron* n = lst.at(randomindex(lst.size()));
						this->AddSynapse(n, b, initWeight, true);
						synsAdded++;

						continue;
					}
					
				}
			}

		}
	}

	printf("Added/removed %d synapses\n", synsAdded);
}






// Perform various tests

#define TEST_CASE( a )  {  int _res = (a); if (!_res) cout << " Failed: " << #a << endl; else cout << "Success: "<< #a << endl; } 

void LANetwork::RunTests()
{

	LANetwork net; 


	uint sz = 20;

	net.CreateNeurons(sz, 10, 'P', &net.pyr_list, -1, 0);
	TEST_CASE( net.pyr_list.size() == sz );

	net.CreateNeurons(sz, 10, 'V', &net.in_pv, -1, 0);
	TEST_CASE( net.in_pv.size() == sz );

	net.CreateNeurons(sz, 10, 'S', &net.bg_list, -1, 0);
	TEST_CASE(net.bg_list.size() == sz);

	for (nrn_iter na = net.pyr_list.begin(); na != net.pyr_list.end(); ++na)
	{
		LANeuron* nrn  = *na;
		TEST_CASE(nrn->type == 'P');
	}

	TEST_CASE(net.synapses.size() ==0 );

	net.ConnectNeurons(net.pyr_list, net.in_pv, false, (float)10.,  1000, 1, 1.0, false);
	TEST_CASE(net.synapses.size()  == 1000);

	net.ConnectNeurons(net.in_pv, net.pyr_list, false, (float)10.,  1000, 1, 1.0, false);
	TEST_CASE(net.synapses.size()  == 2000);

	net.ConnectNeurons(net.bg_list, net.pyr_list, false, (float)10.,  1000, 1, 1.0, false);
	TEST_CASE(net.synapses.size()  == 3000);

	program_input(net.bg_list, 0 , 1000, 100, 1.0, 0);
	for (nrn_iter na = net.bg_list.begin(); na != net.bg_list.end(); ++na)
	{
		LAInput* nrn  = (LAInput*) (*na);
		TEST_CASE(nrn->type == 'S');
		TEST_CASE(nrn->totalSpikes  == 100);
	}


	program_input(net.bg_list, 0 , 1000, 100, 0.0, 0);
	for (nrn_iter na = net.bg_list.begin(); na != net.bg_list.end(); ++na)
	{
		LAInput* nrn  = (LAInput*) (*na);
		TEST_CASE(nrn->totalSpikes  == 100);
	}


	//LAInput* inpA = (LAInput*) net.bg_list[0];
	//LAInput* inpB = (LAInput*) net.bg_list[1];


	net.CalculateDistances();

	for (nrn_iter na = net.neurons.begin(); na != net.neurons.end(); ++na)
		for (nrn_iter nb = net.neurons.begin(); nb != net.neurons.end(); ++nb)
		{
			LANeuron* a = *na, *b=*nb;
			TEST_CASE( (net.distances[pair<int,int>(a->nid,b->nid)] < 1.415) );
		}


}









