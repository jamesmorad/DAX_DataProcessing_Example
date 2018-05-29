#include <stdlib.h>
#include <iostream>
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>

using namespace std;

void S1Filter();
void S2Filter(int inputArr[],int len,Double_t S2filteredArr[]);
void findS2(int inputArr[], Double_t baseline, int len, vector<Int_t>& start,vector<Int_t>& end, vector<Double_t>& area, vector<Double_t>& height, vector<Int_t>& peak,vector<Int_t>& right,vector<Int_t>& left,vector<Double_t>& statWidth,Double_t* S2FilterArr, int event);
void findS1(int inputArr[],Double_t baseline, int len, vector<Int_t>& start,vector<Int_t>& end, vector<Double_t>& area, vector<Double_t>& height, vector<Int_t>& peak,vector<Int_t>& right,vector<Int_t>& left,vector<Double_t>& statWidth,vector<Int_t>& NOTstart,vector<Int_t>& NOTend);

int g_entry = 0;
int g_event = 0;


class RQs
{
	public:
		RQs(bool PIN, TTree* events);
		~RQs();
		void ZeroOut();
		Double_t baselineAvgStart;
		Double_t baselineAvgEnd;
		Double_t baselineStdv;
		Int_t pulseSaturate;
		Double_t minimumVoltage;
		Double_t maximumVoltage;

		vector<Double_t> S1_area_mVns;
		vector<Double_t> S1_height_mV;

		vector<Double_t> S2_area_mVns;
		vector<Double_t> S2_height_mV;

		vector<Int_t> S1_tstart_ns;
		vector<Int_t> S1_t50r_ns;
		vector<Int_t> S1_tpeak_ns;
		vector<Int_t> S1_t50l_ns;
		vector<Int_t> S1_tend_ns;
		vector<Double_t> S1_stat_width;


		vector<Int_t> S2_tstart_ns;
		vector<Int_t> S2_t50r_ns;
		vector<Int_t> S2_tpeak_ns;
		vector<Int_t> S2_t50l_ns;
		vector<Int_t> S2_tend_ns;
		vector<Double_t> S2_stat_width;

		Int_t driftTime_ns;
		Int_t eventNum;
		Int_t numS1;
		Int_t numS2;

		Double_t S2FilterArr[8192];

};
RQs::RQs(bool PIN,TTree* events){

	if(PIN){
		//PINevents stuff
		events->Branch("baselineAvgStart_mV",&this->baselineAvgStart);
		events->Branch("baselineAvgEnd_mV",&this->baselineAvgEnd);
		events->Branch("baselineStdv_mV",&this->baselineStdv);
		events->Branch("maximumVoltage_mV",&this->maximumVoltage);
		events->Branch("minimumVoltage_mV",&this->minimumVoltage);
		events->Branch("pulseSaturate_0NO",&this->pulseSaturate);
	}
	else{
		events->Branch("baselineAvgStart_mV",&this->baselineAvgStart);
		events->Branch("baselineAvgEnd_mV",&this->baselineAvgEnd);
		events->Branch("baselineStdv_mV",&this->baselineStdv);
		events->Branch("pulseSaturate_0NO",&this->pulseSaturate);
		events->Branch("minimumVoltage_mV",&this->minimumVoltage);
		events->Branch("maximumVoltage_mV",&this->maximumVoltage);

		events->Branch("S1_area_mVns",&this->S1_area_mVns);
		events->Branch("S1_height_mV",&this->S1_height_mV);

		events->Branch("S2_area_mVns",&this->S2_area_mVns);
		events->Branch("S2_height_mV",&this->S2_height_mV);


		events->Branch("S1_tstart_ns",&this->S1_tstart_ns);
		events->Branch("S1_t50r_ns",&this->S1_t50r_ns);
		events->Branch("S1_tpeak_ns",&this->S1_tpeak_ns);
		events->Branch("S1_t50l_ns",&this->S1_t50l_ns);
		events->Branch("S1_tend_ns",&this->S1_tend_ns);
		events->Branch("S1_stat_width",&this->S1_stat_width);
		

		events->Branch("S2_tstart_ns",&this->S2_tstart_ns);
		events->Branch("S2_t50r_ns",&this->S2_t50r_ns);
		events->Branch("S2_tpeak_ns",&this->S2_tpeak_ns);
		events->Branch("S2_t50l_ns",&this->S2_t50l_ns);
		events->Branch("S2_tend_ns",&this->S2_tend_ns);	
		events->Branch("eventNum",&this->eventNum);
		events->Branch("S2_stat_width",&this->S2_stat_width);	

		events->Branch("driftTime_ns",&this->driftTime_ns);

		events->Branch("numS1",&this->numS1);
		events->Branch("numS2",&this->numS2);
	}


};
RQs::~RQs(){};
void RQs::ZeroOut(){

	this->baselineAvgStart=0;
	this->baselineAvgEnd=0;
	this->pulseSaturate=0;
	this->minimumVoltage=1e4;
	this->maximumVoltage=-1e4;
	this->baselineStdv = 0;
	this->driftTime_ns = -1; //ADDED (otherwise the previous driftTime_ns gets used again)


}



int main(int argc, char *argv[]){


	if (argc < 3){
		cout << "Pass a ROOT file to process and output filename" << endl;
		return 1;
	}


	////////////////////////////////////////
	////////////////////////////////////////
	////////////////////////////////////////
	//
	//
	//BASIC IO SECTION

	int nSamples;
	int nChannels;
	TString buff;
	int** rawChannelWaveforms = NULL;

	TString inputFilename = argv[1];
	TString outputFilename = argv[2];


	//Open output file and set branches
	auto outputDataFile = new TFile(outputFilename,"RECREATE");
	auto outputEventsLG =  new TTree("eventsLG","eventsLG");
	auto outputEventsHG =  new TTree("eventsHG","eventsHG");
	auto pinEvents = new TTree("PINevents","PINevents");

	auto inputDataFile = new TFile(inputFilename);
	auto events = (TTree*)inputDataFile->Get("events");
	auto runInfo = (TTree*)inputDataFile->Get("runInfo");
	runInfo->SetBranchAddress("nSamples",&nSamples);
	runInfo->SetBranchAddress("nChannels",&nChannels);		
	//Get entry so we can load number of samples and channels into this class
	runInfo->GetEntry(0);
	//Create number of channels using the 2D array/pointer business
	rawChannelWaveforms = new int*[nChannels];
	for(int channel=0;channel<nChannels;channel++){
		//create correct length for each channel array based on number of samples in data file
		//then set the branch address for each channel
		rawChannelWaveforms[channel] = new int[nSamples];
		buff.Form("adcCounts%d",channel);
		events->SetBranchAddress(buff,rawChannelWaveforms[channel]);
	}
	events->SetBranchStatus("*",0);
	
		events->SetBranchStatus("adcCounts0",1);
		events->SetBranchStatus("adcCounts1",1);
	
	if(nChannels>3){
		events->SetBranchStatus("adcCounts4",1);
	}

	
	////////////////////////////////////////
	////////////////////////////////////////
	// ~~ENCAPSULATION~~
	////////////////////////////////////////
	////////////////////////////////////////

	RQs ch0(false,outputEventsHG);
	RQs ch1(false,outputEventsLG);
	RQs ch4(true,pinEvents);

	

	int numEntries = events->GetEntries();
	cout << "There are " << numEntries << " entries. nSamples is " << nSamples << endl;

	for(int entry = 0; entry < numEntries; entry++){

		if(entry % 5000 == 0){ cout << "Processing entry " << entry << " out of " << numEntries << endl;}
		events->GetEntry(entry);
		ch0.ZeroOut();
		ch1.ZeroOut();
		if(nChannels>3){ch4.ZeroOut();}

		
		for(int sample=0;sample<25;sample++){
			ch0.baselineAvgStart += rawChannelWaveforms[0][sample]/25.;
			ch1.baselineAvgStart += rawChannelWaveforms[1][sample]/25.;
			if(nChannels>3){ch4.baselineAvgStart += rawChannelWaveforms[4][sample]/25.;}
		}
    
		//baselineAvgStart*=.122; //CHANGED baseline is used in ADC counts for findS2 findS1

		for(int sample=nSamples-25;sample<nSamples;sample++){
			ch0.baselineAvgEnd += rawChannelWaveforms[0][sample]/25.;
			ch1.baselineAvgEnd += rawChannelWaveforms[1][sample]/25.;
			if(nChannels>3){ch4.baselineAvgEnd += rawChannelWaveforms[4][sample]/25.;}
		}

		ch0.baselineAvgEnd*=.122;
		ch1.baselineAvgEnd*=.122;
		if(nChannels>3){ch4.baselineAvgEnd*=.122;}
		

		for(int sample=0;sample<25;sample++){
			ch0.baselineStdv += (rawChannelWaveforms[0][sample]-ch0.baselineAvgStart)*(rawChannelWaveforms[0][sample]-ch0.baselineAvgStart)/25.;
			ch1.baselineStdv += (rawChannelWaveforms[1][sample]-ch1.baselineAvgStart)*(rawChannelWaveforms[1][sample]-ch1.baselineAvgStart)/25.;
			if(nChannels>3){ch4.baselineStdv += (rawChannelWaveforms[4][sample]-ch4.baselineAvgStart)*(rawChannelWaveforms[4][sample]-ch4.baselineAvgStart)/25.;}
		}

		ch0.baselineStdv = TMath::Sqrt(ch0.baselineStdv);
		ch0.baselineStdv*=.122;
		ch1.baselineStdv = TMath::Sqrt(ch1.baselineStdv);
		ch1.baselineStdv*=.122;
		if(nChannels>3){ch4.baselineStdv = TMath::Sqrt(ch4.baselineStdv);
		ch4.baselineStdv*=.122;}
		
		for(int sample=0;sample<nSamples;sample++){

			////////////////////////////////////////
			//Find minimum voltage
			////////////////////////////////////////
			if(rawChannelWaveforms[0][sample] < ch0.minimumVoltage){
				ch0.minimumVoltage = rawChannelWaveforms[0][sample];
			}
			if(rawChannelWaveforms[1][sample] < ch1.minimumVoltage){
				ch1.minimumVoltage = rawChannelWaveforms[1][sample];
			}
			if(nChannels>3){if(rawChannelWaveforms[4][sample] < ch4.minimumVoltage){
				ch4.minimumVoltage = rawChannelWaveforms[4][sample];
			}}
			////////////////////////////////////////
			//Find maximum voltage
			////////////////////////////////////////
			if(rawChannelWaveforms[0][sample] > ch0.maximumVoltage){
				ch0.maximumVoltage = rawChannelWaveforms[0][sample];
			}
			if(rawChannelWaveforms[1][sample] > ch1.maximumVoltage){
				ch1.maximumVoltage = rawChannelWaveforms[1][sample];
			}
			if(nChannels>3){if(rawChannelWaveforms[4][sample] > ch4.maximumVoltage){
				ch4.maximumVoltage = rawChannelWaveforms[4][sample];
			}}
			////////////////////////////////////////
			//Saturation flag
			////////////////////////////////////////
			if(rawChannelWaveforms[0][sample] > 8191){
				ch0.pulseSaturate = 1;
			}
			if(rawChannelWaveforms[1][sample] > 8191){
				ch1.pulseSaturate = 1;
			}
			if(nChannels>3){if(rawChannelWaveforms[4][sample] > 8191){
				ch4.pulseSaturate = 1;
			}}

			ch0.S2FilterArr[sample]= 0;
			ch1.S2FilterArr[sample]= 0;
			if(nChannels>3){ch4.S2FilterArr[sample]= 0;}
		}

		ch0.minimumVoltage*=.122;
		ch1.minimumVoltage*=.122;
		if(nChannels>3){ch4.minimumVoltage*=.122;}
		
		ch0.maximumVoltage*=.122;
		ch1.maximumVoltage*=.122;
		if(nChannels>3){ch4.maximumVoltage*=.122;}

		
		findS2(rawChannelWaveforms[0], ch0.baselineAvgStart, nSamples, ch0.S2_tstart_ns,ch0.S2_tend_ns, ch0.S2_area_mVns, ch0.S2_height_mV,ch0.S2_tpeak_ns,ch0.S2_t50r_ns,ch0.S2_t50l_ns,ch0.S2_stat_width,ch0.S2FilterArr, entry);
		findS2(rawChannelWaveforms[1], ch1.baselineAvgStart, nSamples, ch1.S2_tstart_ns,ch1.S2_tend_ns, ch1.S2_area_mVns, ch1.S2_height_mV,ch1.S2_tpeak_ns,ch1.S2_t50r_ns,ch1.S2_t50l_ns,ch1.S2_stat_width,ch1.S2FilterArr, entry);
		
		findS1(rawChannelWaveforms[0],ch0.baselineAvgStart, nSamples, ch0.S1_tstart_ns,ch0.S1_tend_ns, ch0.S1_area_mVns, ch0.S1_height_mV,ch0.S1_tpeak_ns,ch0.S1_t50r_ns,ch0.S1_t50l_ns,ch0.S1_stat_width,ch0.S2_tstart_ns,ch0.S2_tend_ns);
		findS1(rawChannelWaveforms[1],ch1.baselineAvgStart, nSamples, ch1.S1_tstart_ns,ch1.S1_tend_ns, ch1.S1_area_mVns, ch1.S1_height_mV,ch1.S1_tpeak_ns,ch1.S1_t50r_ns,ch1.S1_t50l_ns,ch1.S1_stat_width,ch1.S2_tstart_ns,ch1.S2_tend_ns);

		ch0.numS2 = ch0.S2_tstart_ns.size();
		ch1.numS2 = ch1.S2_tstart_ns.size();
		ch0.numS1 = ch0.S1_tstart_ns.size();
		ch1.numS1 = ch1.S1_tstart_ns.size();
		ch0.eventNum = g_event;
		ch1.eventNum = g_event;
		g_event++;

		if(ch0.numS1==1 && ch0.numS2==1){
			ch0.driftTime_ns = ch0.S2_tstart_ns.at(0) - ch0.S1_tstart_ns.at(0);
		}
		if(ch1.numS1==1 && ch1.numS2==1){
			ch1.driftTime_ns = ch1.S2_tstart_ns.at(0) - ch1.S1_tstart_ns.at(0);
		}

		outputEventsLG->Fill();
		outputEventsHG->Fill();
		if(nChannels>3){pinEvents->Fill();}
		
	}


	outputDataFile->cd();
	outputEventsLG->Write();
	outputEventsHG->Write();
	if(nChannels>3){pinEvents->Write();}
	outputDataFile->Close();

	inputDataFile->Close();



	return 0;
	
}


void S2Filter(int inputArr[], int len, Double_t* S2filteredArr){

	Double_t s1windowWidth = 50; //10 samples 100 ns
	Double_t s2windowWidth = 200; //2 us

	Double_t s1Half = (s1windowWidth-1)/2.;
	Double_t s2Half = (s2windowWidth-1)/2.;

	Double_t S1_maxAreaVal;
	
	Double_t outputArrS1[len];
	Double_t outputArrS2[len];

	for(int i=0;i<len;i++){
		outputArrS1[i] =0;
		outputArrS2[i] =0;
		//s1
		for(int j=max(Double_t(0),i-s1Half);j<min(Double_t(len),i+s1Half+1);j++){
			outputArrS1[i] += inputArr[j];
		}
		//s2
		for(int j=max(Double_t(0),i-s2Half);j<min(Double_t(len),i+s2Half+1);j++){
			outputArrS2[i] += inputArr[j];
		}
	}

	for(int i=0;i<len;i++)
  {
		S1_maxAreaVal = 0.;
		for(int j=max(Double_t(0),(i-s2Half+s1Half)); j < min(Double_t(len), (i+s2Half-s1Half+1)); j++)
    {
			if (outputArrS1[j]>S1_maxAreaVal)
				S1_maxAreaVal = outputArrS1[j];
		}

		S2filteredArr[i] = max(Double_t(0.),outputArrS2[i]-S1_maxAreaVal);
	}	
}




bool CheckStartBounds(int i, vector<Int_t>& NOTstart,vector<Int_t>& NOTend){

	bool flag = true;
	for(unsigned int j=0; j<NOTstart.size();j++){
		if ( i*10 >= NOTstart.at(j) &&  i*10 <= NOTend.at(j)){
			flag = false;
		}
	}

	return flag;

}

void findS1(int inputArr[],Double_t baseline, int len, vector<Int_t>& start,vector<Int_t>& end, vector<Double_t>& area, vector<Double_t>& height, vector<Int_t>& peak,vector<Int_t>& right,vector<Int_t>& left,vector<Double_t>& statWidth,vector<Int_t>& NOTstart,vector<Int_t>& NOTend){
	start.clear();
	end.clear();
	area.clear();
	height.clear();
	peak.clear();
	right.clear();
	left.clear();
	statWidth.clear();


	bool foundPulse=false;
	Int_t start_t=0;
	Int_t end_t;
	Double_t area_t=0;
	Double_t height_t=0;
	Int_t peak_t=0;
	Int_t t50r =0;
	Int_t t50l =0;
	Double_t statWidth_t=0;
	Double_t sigma_t2=0;
	Double_t sigma_t=0;

	Double_t moving_area_calculator=0;
	bool found50r;
	bool found50l;

	for(int i=0;i<len;i++){

		
		if( CheckStartBounds(i,NOTstart,NOTend) ) {

			if( inputArr[i] - baseline > 40 && !foundPulse){ //CHANGE: 40 ADCC based on sample to baseline plot
				foundPulse=true;
				start_t = max(0,i-4); //changed 3 to 4 to accomadate 
				
				for(int k=max(0,i-4);k<i;k++){
					area_t += (inputArr[k]- baseline)*10*.122;
				}
				area_t += (inputArr[i]- baseline)*10*.122;
				if( (inputArr[i]- baseline)*.122 > height_t){
					height_t = (inputArr[i]- baseline)*.122;
					peak_t=i;
				}
			}

			if(foundPulse){
				area_t += (inputArr[i]- baseline)*10*.122;
				if( (inputArr[i]- baseline)*.122 > height_t){
					height_t = (inputArr[i]- baseline)*.122;
					peak_t=i;
				}
			}
			if(foundPulse && inputArr[i] - baseline < 10){ 

				end_t = min(i+3,len); 
				for(int k=i+1;k<end_t;k++){
					area_t += (inputArr[k]- baseline)*10*.122;
				}
				moving_area_calculator=0;
				found50r=false;
				found50l=false;


				for(int k=start_t;k<=end_t;k++){

					sigma_t2 += (k-start_t)*(k-start_t)*(inputArr[k]- baseline)/(end_t-start_t);
					sigma_t += (k-start_t)*(inputArr[k]- baseline)/(end_t-start_t);

					moving_area_calculator += (inputArr[k]- baseline)*10*.122;
					if(moving_area_calculator > 0.25*area_t && !found50r){
						t50r = k;
						found50r=true;
					}
					if(moving_area_calculator > 0.75*area_t && !found50l){
						t50l = k;
						found50l=true;
					}
				}

				sigma_t2 = sigma_t2/area_t;
				sigma_t = sigma_t/area_t;

				statWidth_t = TMath::Sqrt(sigma_t2 - sigma_t*sigma_t);

        //cout << "\tS1 start: " << start_t << " end: " << end_t << endl; 

				area.push_back(area_t);
				start.push_back(10*start_t);
				end.push_back(10*end_t);
				height.push_back(height_t);
				peak.push_back(10*peak_t);
				right.push_back(10*t50r);
				left.push_back(10*t50l);
				statWidth.push_back(10*statWidth_t);

				area_t=0;
				start_t=0;
				end_t=0;
				height_t=0;
				peak_t=0;
				t50r=0;
				t50l=0;
				sigma_t2=0;
				sigma_t=0;
				statWidth_t=0;
				foundPulse = false;
			}
		}
	}
}

void findS2(int inputArr[],Double_t baseline, int len, vector<Int_t>& start,vector<Int_t>& end, vector<Double_t>& area, vector<Double_t>& height, vector<Int_t>& peak,vector<Int_t>& right,vector<Int_t>& left,vector<Double_t>& statWidth,Double_t* S2FilterArr, int event){

	start.clear();
	end.clear();
	area.clear();
	height.clear();
	peak.clear();
	right.clear();
	left.clear();
	statWidth.clear();
	//Double_t S2filteredArr[len];
	S2Filter(inputArr,len,S2FilterArr);

	bool foundPulse=false;
	Int_t start_t=0;
	Int_t end_t=0;
	Double_t area_t=0;
	Double_t height_t=0;


	Int_t peak_t=0;
	Int_t t50r =0;
	Int_t t50l =0;
	Double_t statWidth_t=0;
	Double_t sigma_t2=0;
	Double_t sigma_t=0;

	Double_t moving_area_calculator=0;
	bool found50r;
	bool found50l;

	for(int i=0; i<len;i++){

		if( S2FilterArr[i] > 2500. && !foundPulse){
			foundPulse=true;
			//adding in a 100ns pretrigger for the S2
			start_t = max(i-10,0);
			for(int pretrigger = start_t;pretrigger<=i;pretrigger++){
				area_t += (inputArr[pretrigger]- baseline)*10*.122;
			}

			if( (inputArr[i]- baseline)*.122 > height_t){
				height_t = (inputArr[i]- baseline)*.122;
				peak_t = i;
			}
		}

		if (foundPulse){
			area_t += (inputArr[i]- baseline)*10*.122;
			if( (inputArr[i]- baseline)*.122 > height_t){
				height_t = (inputArr[i]- baseline)*.122;
				peak_t = i;
			}
		}

		if(foundPulse && S2FilterArr[i] < 2000.){
			end_t = min(i+50,len);

			for(int j=i+1; j < end_t; j++){
				area_t += (inputArr[j]- baseline)*10*.122;
			}
			moving_area_calculator=0;
			found50r=false;
			found50l=false;

			for(int k=start_t;k<=end_t;k++){

				sigma_t2 += (k-start_t)*(k-start_t)*(inputArr[k]- baseline)/(end_t-start_t);
				sigma_t += (k-start_t)*(inputArr[k]- baseline)/(end_t-start_t);

				moving_area_calculator += (inputArr[k]- baseline)*10*.122;
				if(moving_area_calculator > 0.25*area_t && !found50r){
					t50r = k;
					found50r=true;
				}
				if(moving_area_calculator > 0.75*area_t && !found50l){
					t50l = k;
					found50l=true;
				}
			}

			sigma_t2 = sigma_t2/area_t;
			sigma_t = sigma_t/area_t;

			statWidth_t = TMath::Sqrt(sigma_t2 - sigma_t*sigma_t);

      //cout << "entry " << g_entry << " (event " << event << ") has S2 area " << area_t << (area_t < 20e3? "  ***" : "") << endl;;
      //cout << "\tS2 start: " << start_t << " end: " << end_t << endl;

      g_entry++;

			area.push_back(area_t);
			start.push_back(10*start_t);
			end.push_back(10*end_t);
			height.push_back(height_t);
			peak.push_back(10*peak_t);
			right.push_back(10*t50r);
			left.push_back(10*t50l);
			statWidth.push_back(10*statWidth_t);

			area_t=0;
			start_t=0;
			end_t=0;
			height_t=0;
			peak_t=0;
			t50r=0;
			t50l=0;
			sigma_t2=0;
			sigma_t=0;
			statWidth_t=0;
			foundPulse = false;
		}
	}
}


