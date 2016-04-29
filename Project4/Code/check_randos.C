#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <map>

#include "randos.C"
#include "classes.C"

using namespace std;

void freq(string s, int* pt){
  /*
    Return a vector of the number of times each digit appears in the sequence
  */

  for(std::string::iterator it=s.begin(); it!=s.end(); ++it){
    stringstream temp1;
    temp1 << *it;
    string n;
    temp1 >> n;
    pt[atoi(n.c_str())]++;
  }

}

void serial(string s, int* pt){
  /*
    Returns a vector of the number of times each 2-digit pair appears in the sequence
  */

  for(unsigned int i=0;i<s.size();i++){
    stringstream temp;
    if(i==s.size()-1)
      temp<<s.at(i)<<s.at(0);
    else
      temp<<s.at(i)<<s.at(i+1);
    string hold;
    temp >> hold;
    int n = atoi(hold.c_str());
    pt[n]++;
  }

}

void poker(string s, int *pt){
  /*
    Returns a vector of the number of times each poker hand occurs (5, 4, 3, 2, full, two pair, other)
  */

  for(int it=0;it<s.size();it++){
    
    string hand;
    stringstream temp;
    temp<<s.at(it)<<s.at((it+1)%s.size())<<s.at((it+2)%s.size())<<s.at((it+3)%s.size())<<s.at((it+4)%s.size());
    temp>>hand;
    map<char,int> ct;
    vector<char> imports;
    for(std::string::iterator i=hand.begin(); i!=hand.end(); i++){
      char hold = *i;
      if(ct.find(hold)==ct.end()){
	ct[hold]=1;
	imports.push_back(hold);
      }
      else
	ct[hold]++;
    }
    if(ct.size()==5)
      pt[6]++;
    else if(ct.size()==1)
      pt[0]++;
    else if(ct.size()==2){
      //4 of a kind or full house
      if(ct[imports.at(0)]==4 || ct[imports.at(0)]==1)
	pt[1]++;
      else
	pt[4]++;
    }
    else if(ct.size()==4)
      //One pair
      pt[3]++;
    else{
      //3 of a kind or two pair
      if(ct[imports.at(0)]==3 || ct[imports.at(1)]==3 || ct[imports.at(2)]==3)
	pt[2]++;
      else
	pt[5]++;
    }
  }

}

void gap(string s, int *pt){
  /*
    Returns a vector of the number of times each gap size occurs between 0's (0,1,...,15,16-20,21-25,26+)
  */

  int ct = -1;
  for(int it=0;it<s.size();it++){
    
    if(s.at(it)=='0'){
      if(ct==-1)
	ct+=1;
      else if(ct < 16){
	pt[ct]++; ct=0;}
      else if(ct < 21){
	pt[16]++; ct=0;}
      else if(ct < 26){
	pt[17]++; ct=0;}
      else{
	pt[18]++; ct=0;}
    }
    else if(ct!=-1)
      ct++;
  }

}

double chi2(int *obs, double *exp, int sz){
  /*
    Returns a double of chi^2
  */

  double to_ret = 0;

  for(int i=0;i<sz;i++){
    //cout<<"obs: "<<obs[i]<<"; exp: "<<exp[i]<<endl;
    if(obs[i]!=0)
      to_ret+=1.0*pow((obs[i]-exp[i]),2)/obs[i];
    else
      cout<<"ERROR: "<<i<<endl;
  }

  return to_ret;
}

double mean(string s){
  /*
    Returns a double of the mean value of the digits in the string
  */

  int ct=0;
  int sum = 0;
  for(std::string::iterator it=s.begin(); it!=s.end(); it++){
    stringstream temp;
    temp << *it;
    string n;
    temp >> n;
    sum+=atoi(n.c_str());
    ct++;
  }

  return 1.0*sum/ct;

}

double std_dev(string s, double mu){
  /*
    Returns a double of the standard deviation from the mean of the digits in the string
  */

  int ct=0;
  int sum = 0;
  for(std::string::iterator it=s.begin(); it!=s.end(); it++){
    stringstream temp;
    temp << *it;
    string n;
    temp >> n;
    sum+=pow(atoi(n.c_str()),2);
    ct++;
  }

  double mean2 = 1.0*sum/ct;

  return sqrt(abs(mean2-pow(mu,2)));

}

void check_rando0(double len){
  /*
    Generates a sequence of digits using ran0 and runs various tests.
  */

  ofstream myfile;
  myfile.open("benchmarks/ran0_check.txt"); 

  long int seed = 1000;
  double ran = ran0(&seed);
  //cout<<ran<<endl;

  string s;

  while(s.size()<len){
    string temp = to_string(ran);
    /*if(temp.at(0)=='0'){
      stringstream hold;
      for(std::string::iterator it=temp.begin()+2; it!=temp.end(); it++){
	char n = *it;
	hold<<n;
      }
      hold>>temp;
      for(std::string::iterator it=temp.begin(); it!=temp.end(); it++){
	stringstream temp1;
	temp1 << *it;
	string n;
	temp1 >> n;
      }
    }
    while(temp.size() < 6)
      temp+="0";
    for(int i=0;i<temp.size();i++){
      if(s.size()==len)
	continue;
      s+=temp.at(i);
      }*/
    if(temp.at(0)=='0')
      s+=temp.at(2);
    else
      s+=temp.at(0);
    ran = ran0(&seed);
    //cout<<ran<<endl;
  }

  //Run frequency test
  myfile<<endl<<"Ran0"<<endl<<endl;
  myfile<<"Frequency"<<endl;
  int *frequency_test;
  frequency_test = new int[10];
  for(int i=0;i<10;i++)
    frequency_test[i]=0;

  freq(s,frequency_test);

  double *frequency_exp;
  frequency_exp = new double[10];
  for(int i=0;i<10;i++)
    frequency_exp[i]=len*0.1;

  double chi_frequency = chi2(frequency_test,frequency_exp,10);

  //Run serial test
  myfile<<"serial"<<endl;
  int *serial_test;
  serial_test = new int[100];
  for(int i=0;i<100;i++)
    serial_test[i]=0;

  serial(s,serial_test);

  double *serial_exp;
  serial_exp = new double[100];
  for(int i=0;i<100;i++)
    serial_exp[i]=len*0.01;

  double chi_serial = chi2(serial_test,serial_exp,100);

  //Run poker test
  myfile<<"poker"<<endl;
  int *poker_test;
  poker_test = new int[7];
  for(int i=0;i<7;i++)
    poker_test[i]=0;

  poker(s,poker_test);

  double *poker_exp;
  poker_exp = new double[7];
  poker_exp[0] = 0.0001*len; poker_exp[1] = 0.0045*len; poker_exp[2] = 0.072*len;
  poker_exp[3] = 0.504*len; poker_exp[4] = 0.009*len; poker_exp[5] = 0.108*len; poker_exp[6] = 0.3024*len;

  double chi_poker = chi2(poker_test,poker_exp,7);

  //Run gap test
  myfile<<"gap"<<endl;
  int *gap_test;
  gap_test = new int[19];
  for(int i=0;i<19;i++)
    gap_test[i]=0;

  gap(s,gap_test);

  double *gap_exp;
  gap_exp = new double[19];
  gap_exp[0]=len*48.3/5000; gap_exp[1]=len*43.5/5000; gap_exp[2]=len*39.1/5000; gap_exp[3]=len*35.2/5000; 
  gap_exp[4]=len*31.7/5000;gap_exp[5]=len*28.5/5000; gap_exp[6]=len*25.7/5000; gap_exp[7]=len*23.1/5000;
  gap_exp[8]=len*20.8/5000; gap_exp[9]=len*18.7/5000; gap_exp[10]=len*16.8/5000; gap_exp[11]=len*15.2/5000;
  gap_exp[12]=len*13.6/5000; gap_exp[13]=len*12.3/5000; gap_exp[14]=len*11.0/5000; gap_exp[15]=len*9.9/5000; 
  gap_exp[16]=len*36.7/5000; gap_exp[17]=len*21.6/5000; gap_exp[18]=len*31.2/5000;

  double chi_gap = chi2(gap_test,gap_exp,19);

  //Compute mean
  double mu = mean(s);

  //Compute std dev
  double sigma = std_dev(s,mu);

  myfile<<endl<<"Ran0: "<<endl<<endl;//s<<endl<<endl;
  myfile<<"Mean: "<<mu<<endl;
  myfile<<"Sigm: "<<sigma<<endl;
  myfile<<"Freq chi: "<<chi_frequency<<endl;
  myfile<<"Seri chi: "<<chi_serial<<endl;
  myfile<<"Poke chi: "<<chi_poker<<endl;
  myfile<<"Gap  chi: "<<chi_gap<<endl<<endl;

  delete [] frequency_test;
  delete [] frequency_exp;
  delete [] serial_test;
  delete [] serial_exp;
  delete [] poker_test;
  delete [] poker_exp;
  delete [] gap_test;
  delete [] gap_exp;

  myfile.close();

}

void check_rando1(double len){
  /*
    Generates a sequence of digits using ran1 and runs various tests.
  */

  ofstream myfile;
  myfile.open("benchmarks/ran1_check.txt"); 
  long int seed = 1000;
  double ran = ran1(&seed);
  //cout<<ran<<endl;

  string s;

  while(s.size()<len){
    string temp = to_string(ran);
    /*if(temp.at(0)=='0'){
      /* stringstream hold;
      for(std::string::iterator it=temp.begin()+2; it!=temp.end(); it++){
	char n = *it;
	hold<<n;
      }
      hold>>temp;
      for(std::string::iterator it=temp.begin(); it!=temp.end(); it++){
	stringstream temp1;
	temp1 << *it;
	string n;
	temp1 >> n;
      }
    }
    while(temp.size()<6)
      temp+="0";
    for(int i=0;i<temp.size();i++){
      if(s.size()==len)
	continue;
      s+=temp.at(i);
      }*/
    if(temp.at(0)=='0')
      s+=temp.at(2);
    else
      s+=temp.at(0);

    ran = ran1(&seed);
    //cout<<ran<<endl;
  }

  //Run frequency test
  myfile<<endl<<"Ran1"<<endl<<endl;
  myfile<<"Frequency"<<endl;
  int *frequency_test = new int[10];
  for(int i=0;i<10;i++)
    frequency_test[i]=0;

  freq(s,frequency_test);

  double *frequency_exp = new double[10];
  for(int i=0;i<10;i++)
    frequency_exp[i]=len*0.1;

  double chi_frequency = chi2(frequency_test,frequency_exp,10);

  //Run serial test
  myfile<<"serial"<<endl;
  int *serial_test = new int[100];
  for(int i=0;i<100;i++)
    serial_test[i]=0;

  serial(s,serial_test);

  double *serial_exp = new double[100];
  for(int i=0;i<100;i++)
    serial_exp[i]=len*0.01;

  double chi_serial = chi2(serial_test,serial_exp,100);

  //Run poker test
  myfile<<"poker"<<endl;
  int *poker_test = new int[7];
  for(int i=0;i<7;i++)
    poker_test[i]=0;

  poker(s,poker_test);

  double *poker_exp = new double[7];
  poker_exp[0] = 0.0001*len; poker_exp[1] = 0.0045*len; poker_exp[2] = 0.072*len;
  poker_exp[3] = 0.504*len; poker_exp[4] = 0.009*len; poker_exp[5] = 0.108*len; poker_exp[6] = 0.3024*len;

  double chi_poker = chi2(poker_test,poker_exp,7);

  //Run gap test
  myfile<<"gap"<<endl;
  int *gap_test = new int[19];
  for(int i=0;i<19;i++)
    gap_test[i]=0;

  gap(s,gap_test);

  double *gap_exp = new double[19];
  gap_exp[0]=len*48.3/5000; gap_exp[1]=len*43.5/5000; gap_exp[2]=len*39.1/5000; gap_exp[3]=len*35.2/5000; 
  gap_exp[4]=len*31.7/5000;gap_exp[5]=len*28.5/5000; gap_exp[6]=len*25.7/5000; gap_exp[7]=len*23.1/5000;
  gap_exp[8]=len*20.8/5000; gap_exp[9]=len*18.7/5000; gap_exp[10]=len*16.8/5000; gap_exp[11]=len*15.2/5000;
  gap_exp[12]=len*13.6/5000; gap_exp[13]=len*12.3/5000; gap_exp[14]=len*11.0/5000; gap_exp[15]=len*9.9/5000; 
  gap_exp[16]=len*36.7/5000; gap_exp[17]=len*21.6/5000; gap_exp[18]=len*31.2/5000;

  double chi_gap = chi2(gap_test,gap_exp,19);

  //Compute mean
  double mu = mean(s);

  //Compute std dev
  double sigma = std_dev(s,mu);

  myfile<<endl<<"Ran1: "<<endl<<endl;//s<<endl<<endl;
  myfile<<"Mean: "<<mu<<endl;
  myfile<<"Sigm: "<<sigma<<endl;
  myfile<<"Freq chi: "<<chi_frequency<<endl;
  myfile<<"Seri chi: "<<chi_serial<<endl;
  myfile<<"Poke chi: "<<chi_poker<<endl;
  myfile<<"Gap  chi: "<<chi_gap<<endl<<endl;

  delete [] frequency_test;
  delete [] frequency_exp;
  delete [] serial_test;
  delete [] serial_exp;
  delete [] poker_test;
  delete [] poker_exp;
  delete [] gap_test;
  delete [] gap_exp;

  myfile.close();
}

void check_rando2(double len){
  /*
    Generates a sequence of digits using ran2 and runs various tests.
  */

  ofstream myfile;
  myfile.open("benchmarks/ran2_check.txt"); 
  long int seed = 1000;
  double ran = ran2(&seed);
  //cout<<ran<<endl;

  string s;

  while(s.size()<len){
    string temp = to_string(ran);
    /*if(temp.at(0)=='0'){
      stringstream hold;
      for(std::string::iterator it=temp.begin()+2; it!=temp.end(); it++){
	char n = *it;
	hold<<n;
      }
      hold>>temp;
      for(std::string::iterator it=temp.begin(); it!=temp.end(); it++){
	stringstream temp1;
	temp1 << *it;
	string n;
	temp1 >> n;
      }
    }
    while(temp.size()<6)
      temp+="0";
    for(int i=0;i<temp.size();i++){
      if(s.size()==len)
	continue;
      s+=temp.at(i);
      }*/
    if(temp.at(0)=='0')
      s+=temp.at(2);
    else
      s+=temp.at(0);

    ran = ran2(&seed);
    //cout<<ran<<endl;
  }

  //Run frequency test
  myfile<<endl<<"Ran2"<<endl<<endl;
  myfile<<"Frequency"<<endl;
  int *frequency_test = new int[10];
  for(int i=0;i<10;i++)
    frequency_test[i]=0;

  freq(s,frequency_test);

  double *frequency_exp = new double[10];
  for(int i=0;i<10;i++)
    frequency_exp[i]=len*0.1;

  double chi_frequency = chi2(frequency_test,frequency_exp,10);

  //Run serial test
  myfile<<"serial"<<endl;
  int *serial_test = new int[100];
  for(int i=0;i<100;i++)
    serial_test[i]=0;

  serial(s,serial_test);

  double *serial_exp = new double[100];
  for(int i=0;i<100;i++)
    serial_exp[i]=len*0.01;

  double chi_serial = chi2(serial_test,serial_exp,100);

  //Run poker test
  myfile<<"poker"<<endl;
  int *poker_test = new int[7];
  for(int i=0;i<7;i++)
    poker_test[i]=0;

  poker(s,poker_test);

  double *poker_exp = new double[7];
  poker_exp[0] = 0.0001*len; poker_exp[1] = 0.0045*len; poker_exp[2] = 0.072*len;
  poker_exp[3] = 0.504*len; poker_exp[4] = 0.009*len; poker_exp[5] = 0.108*len; poker_exp[6] = 0.3024*len;

  double chi_poker = chi2(poker_test,poker_exp,7);

  //Run gap test
  myfile<<"gap"<<endl;
  int *gap_test = new int[19];
  for(int i=0;i<19;i++)
    gap_test[i]=0;

  gap(s,gap_test);

  double *gap_exp = new double[19];
  gap_exp[0]=len*48.3/5000; gap_exp[1]=len*43.5/5000; gap_exp[2]=len*39.1/5000; gap_exp[3]=len*35.2/5000; 
  gap_exp[4]=len*31.7/5000;gap_exp[5]=len*28.5/5000; gap_exp[6]=len*25.7/5000; gap_exp[7]=len*23.1/5000;
  gap_exp[8]=len*20.8/5000; gap_exp[9]=len*18.7/5000; gap_exp[10]=len*16.8/5000; gap_exp[11]=len*15.2/5000;
  gap_exp[12]=len*13.6/5000; gap_exp[13]=len*12.3/5000; gap_exp[14]=len*11.0/5000; gap_exp[15]=len*9.9/5000; 
  gap_exp[16]=len*36.7/5000; gap_exp[17]=len*21.6/5000; gap_exp[18]=len*31.2/5000;

  double chi_gap = chi2(gap_test,gap_exp,19);

  //Compute mean
  double mu = mean(s);

  //Compute std dev
  double sigma = std_dev(s,mu);

  myfile<<endl<<"Ran2: "<<endl<<endl;//s<<endl<<endl;
  myfile<<"Mean: "<<mu<<endl;
  myfile<<"Sigm: "<<sigma<<endl;
  myfile<<"Freq chi: "<<chi_frequency<<endl;
  myfile<<"Seri chi: "<<chi_serial<<endl;
  myfile<<"Poke chi: "<<chi_poker<<endl;
  myfile<<"Gap  chi: "<<chi_gap<<endl<<endl;

  delete [] frequency_test;
  delete [] frequency_exp;
  delete [] serial_test;
  delete [] serial_exp;
  delete [] poker_test;
  delete [] poker_exp;
  delete [] gap_test;
  delete [] gap_exp;

  myfile.close();
}

void check_rando3(double len){
  /*
    Generates a sequence of digits using ran3 and runs various tests.
  */

  ofstream myfile;
  myfile.open("benchmarks/ran3_check.txt"); 
  long int seed = 1000;
  double ran = ran3(&seed);
  //cout<<ran<<endl;

  string s;

  while(s.size()<len){
    string temp = to_string(ran);
    /*if(temp.at(0)=='0'){
      stringstream hold;
      for(std::string::iterator it=temp.begin()+2; it!=temp.end(); it++){
	char n = *it;
	hold<<n;
      }
      hold>>temp;
      for(std::string::iterator it=temp.begin(); it!=temp.end(); it++){
	stringstream temp1;
	temp1 << *it;
	string n;
	temp1 >> n;
      }
    }
    while(temp.size()<6)
      temp+="0";
    for(int i=0;i<temp.size();i++){
      if(s.size()==len)
	continue;
      s+=temp.at(i);
      }*/
    if(temp.at(0)=='0')
      s+=temp.at(2);
    else
      s+=temp.at(0);

    ran = ran3(&seed);
    //cout<<ran<<endl;
  }

  //Run frequency test
  myfile<<endl<<"Ran3"<<endl<<endl;
  myfile<<"Frequency"<<endl;
  int *frequency_test = new int[10];
  for(int i=0;i<10;i++)
    frequency_test[i]=0;

  freq(s,frequency_test);

  double *frequency_exp = new double[10];
  for(int i=0;i<10;i++)
    frequency_exp[i]=len*0.1;

  double chi_frequency = chi2(frequency_test,frequency_exp,10);

  //Run serial test
  myfile<<"serial"<<endl;
  int *serial_test = new int[100];
  for(int i=0;i<100;i++)
    serial_test[i]=0;

  serial(s,serial_test);

  double *serial_exp = new double[100];
  for(int i=0;i<100;i++)
    serial_exp[i]=len*0.01;

  double chi_serial = chi2(serial_test,serial_exp,100);

  //Run poker test
  myfile<<"poker"<<endl;
  int *poker_test = new int[7];
  for(int i=0;i<7;i++)
    poker_test[i]=0;

  poker(s,poker_test);

  double *poker_exp = new double[7];
  poker_exp[0] = 0.0001*len; poker_exp[1] = 0.0045*len; poker_exp[2] = 0.072*len;
  poker_exp[3] = 0.504*len; poker_exp[4] = 0.009*len; poker_exp[5] = 0.108*len; poker_exp[6] = 0.3024*len;

  double chi_poker = chi2(poker_test,poker_exp,7);

  //Run gap test
  myfile<<"gap"<<endl;
  int *gap_test = new int[19];
  for(int i=0;i<19;i++)
    gap_test[i]=0;

  gap(s,gap_test);

  double *gap_exp = new double[19];
  gap_exp[0]=len*48.3/5000; gap_exp[1]=len*43.5/5000; gap_exp[2]=len*39.1/5000; gap_exp[3]=len*35.2/5000; 
  gap_exp[4]=len*31.7/5000;gap_exp[5]=len*28.5/5000; gap_exp[6]=len*25.7/5000; gap_exp[7]=len*23.1/5000;
  gap_exp[8]=len*20.8/5000; gap_exp[9]=len*18.7/5000; gap_exp[10]=len*16.8/5000; gap_exp[11]=len*15.2/5000;
  gap_exp[12]=len*13.6/5000; gap_exp[13]=len*12.3/5000; gap_exp[14]=len*11.0/5000; gap_exp[15]=len*9.9/5000; 
  gap_exp[16]=len*36.7/5000; gap_exp[17]=len*21.6/5000; gap_exp[18]=len*31.2/5000;

  double chi_gap = chi2(gap_test,gap_exp,19);

  //Compute mean
  double mu = mean(s);

  //Compute std dev
  double sigma = std_dev(s,mu);

  myfile<<endl<<"Ran3: "<<endl<<endl;//s<<endl<<endl;
  myfile<<"Mean: "<<mu<<endl;
  myfile<<"Sigm: "<<sigma<<endl;
  myfile<<"Freq chi: "<<chi_frequency<<endl;
  myfile<<"Seri chi: "<<chi_serial<<endl;
  myfile<<"Poke chi: "<<chi_poker<<endl;
  myfile<<"Gap  chi: "<<chi_gap<<endl<<endl;

  delete [] frequency_test;
  delete [] frequency_exp;
  delete [] serial_test;
  delete [] serial_exp;
  delete [] poker_test;
  delete [] poker_exp;
  delete [] gap_test;
  delete [] gap_exp;

  myfile.close();
}

void check_rando4(double len){
  /*
    Generates a sequence of digits using ran4 (C++ random) and runs various tests.
  */

  ofstream myfile;
  myfile.open("benchmarks/ran4_check.txt"); 
  srand(14);
  double ran = rand()%10;

  string s = to_string(ran);

  while(s.size()<len){
    string temp = to_string(ran);
    s+=temp;
    /*for(int i=0;i<temp.size();i++){
      if(s.size()==len)
	continue;
      s+=temp.at(i);
      }*/
    ran = rand()%10;
    //cout<<ran<<endl;
  }

  //Run frequency test
  myfile<<endl<<"Ran4"<<endl<<endl;
  myfile<<"Frequency"<<endl;
  int *frequency_test = new int[10];
  for(int i=0;i<10;i++)
    frequency_test[i]=0;

  freq(s,frequency_test);

  double *frequency_exp = new double[10];
  for(int i=0;i<10;i++)
    frequency_exp[i]=len*0.1;

  double chi_frequency = chi2(frequency_test,frequency_exp,10);

  //Run serial test
  myfile<<"serial"<<endl;
  int *serial_test = new int[100];
  for(int i=0;i<100;i++)
    serial_test[i]=0;

  serial(s,serial_test);

  double *serial_exp = new double[100];
  for(int i=0;i<100;i++)
    serial_exp[i]=len*0.01;

  double chi_serial = chi2(serial_test,serial_exp,100);

  //Run poker test
  myfile<<"poker"<<endl;
  int *poker_test = new int[7];
  for(int i=0;i<7;i++)
    poker_test[i]=0;

  poker(s,poker_test);

  double *poker_exp = new double[7];
  poker_exp[0] = 0.0001*len; poker_exp[1] = 0.0045*len; poker_exp[2] = 0.072*len;
  poker_exp[3] = 0.504*len; poker_exp[4] = 0.009*len; poker_exp[5] = 0.108*len; poker_exp[6] = 0.3024*len;

  double chi_poker = chi2(poker_test,poker_exp,7);

  //Run gap test
  myfile<<"gap"<<endl;
  int *gap_test = new int[19];
  for(int i=0;i<19;i++)
    gap_test[i]=0;

  gap(s,gap_test);

  double *gap_exp = new double[19];
  gap_exp[0]=len*48.3/5000; gap_exp[1]=len*43.5/5000; gap_exp[2]=len*39.1/5000; gap_exp[3]=len*35.2/5000; 
  gap_exp[4]=len*31.7/5000;gap_exp[5]=len*28.5/5000; gap_exp[6]=len*25.7/5000; gap_exp[7]=len*23.1/5000;
  gap_exp[8]=len*20.8/5000; gap_exp[9]=len*18.7/5000; gap_exp[10]=len*16.8/5000; gap_exp[11]=len*15.2/5000;
  gap_exp[12]=len*13.6/5000; gap_exp[13]=len*12.3/5000; gap_exp[14]=len*11.0/5000; gap_exp[15]=len*9.9/5000; 
  gap_exp[16]=len*36.7/5000; gap_exp[17]=len*21.6/5000; gap_exp[18]=len*31.2/5000;

  double chi_gap = chi2(gap_test,gap_exp,19);

  //Compute mean
  double mu = mean(s);

  //Compute std dev
  double sigma = std_dev(s,mu);

  myfile<<endl<<"Ran4: "<<endl<<endl;//s<<endl<<endl;
  myfile<<"Mean: "<<mu<<endl;
  myfile<<"Sigm: "<<sigma<<endl;
  myfile<<"Freq chi: "<<chi_frequency<<endl;
  myfile<<"Seri chi: "<<chi_serial<<endl;
  myfile<<"Poke chi: "<<chi_poker<<endl;
  myfile<<"Gap  chi: "<<chi_gap<<endl<<endl;

  delete [] frequency_test;
  delete [] frequency_exp;
  delete [] serial_test;
  delete [] serial_exp;
  delete [] poker_test;
  delete [] poker_exp;
  delete [] gap_test;
  delete [] gap_exp;

  myfile.close();
}

void check_randos(){
  /*
    Main function - calls all other functions
  */

  check_rando0(1000000);
  check_rando1(1000000);
  check_rando2(1000000);
  check_rando3(1000000);
  check_rando4(1000000);
}
