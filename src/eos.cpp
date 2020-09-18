#include "eos.h"


void EoS::is_file_exist(fstream& file)
{
  if(!file){cout<<"Files for EoS not found, code terminated. exiting ...."<<endl; exit(1);}
  //else{cout<<"File found ..."<<endl;}
}

EoS::~EoS()
{
    gsl_spline_free (spline_1);
    gsl_spline_free (spline_2);
    gsl_spline_free (spline_3);
    gsl_spline_free (spline_4);
    gsl_spline_free (spline_5);
    gsl_spline_free (spline_6);
    gsl_spline_free (spline_7);

    gsl_spline_free (spline_prs_1);
    gsl_spline_free (spline_prs_2);
    gsl_spline_free (spline_prs_3);
    gsl_spline_free (spline_prs_4);
    gsl_spline_free (spline_prs_5);
    gsl_spline_free (spline_prs_6);
    gsl_spline_free (spline_prs_7);

    gsl_interp_accel_free (acc);

}

EoS::EoS()
{
  
  infile1_d.open("EoS_s95p-v1/s95p-v1_dens1.dat",ios::in);
  infile2_d.open("EoS_s95p-v1/s95p-v1_dens2.dat",ios::in);
  infile3_d.open("EoS_s95p-v1/s95p-v1_dens3.dat",ios::in);
  infile4_d.open("EoS_s95p-v1/s95p-v1_dens4.dat",ios::in);
  infile5_d.open("EoS_s95p-v1/s95p-v1_dens5.dat",ios::in);
  infile6_d.open("EoS_s95p-v1/s95p-v1_dens6.dat",ios::in);
  infile7_d.open("EoS_s95p-v1/s95p-v1_dens7.dat",ios::in);

  infile1_t.open("EoS_s95p-v1/s95p-v1_par1.dat",ios::in);    
  infile2_t.open("EoS_s95p-v1/s95p-v1_par2.dat",ios::in);    
  infile3_t.open("EoS_s95p-v1/s95p-v1_par3.dat",ios::in);    
  infile4_t.open("EoS_s95p-v1/s95p-v1_par4.dat",ios::in);    
  infile5_t.open("EoS_s95p-v1/s95p-v1_par5.dat",ios::in);    
  infile6_t.open("EoS_s95p-v1/s95p-v1_par6.dat",ios::in);    
  infile7_t.open("EoS_s95p-v1/s95p-v1_par7.dat",ios::in);    

  is_file_exist(infile1_d);
  is_file_exist(infile2_d);
  is_file_exist(infile3_d);
  is_file_exist(infile4_d);
  is_file_exist(infile5_d);
  is_file_exist(infile6_d);
  is_file_exist(infile7_d);

  is_file_exist(infile1_t);
  is_file_exist(infile2_t);
  is_file_exist(infile3_t);
  is_file_exist(infile4_t);
  is_file_exist(infile5_t);
  is_file_exist(infile6_t);
  is_file_exist(infile7_t);


  infile1_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> e0_1;
  delete iss;
  
  infile1_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> de_1 >> ne_1 ;
  delete iss;

  infile2_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> e0_2;
  delete iss;
  
  infile2_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> de_2 >> ne_2 ;
  delete iss;

  infile3_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> e0_3;
  delete iss;
  
  infile3_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> de_3 >> ne_3 ;
  delete iss;

  infile4_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> e0_4;
  delete iss;
  
  infile4_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> de_4 >> ne_4 ;
  delete iss;

  infile5_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> e0_5;
  delete iss;
  
  infile5_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> de_5 >> ne_5 ;
  delete iss;

  infile6_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> e0_6;
  delete iss;
  
  infile6_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> de_6 >> ne_6 ;
  delete iss;

  infile7_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> e0_7;
  delete iss;
  
  infile7_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> de_7 >> ne_7 ;
  delete iss;


  
  for(int i=ne_1-1; i>=0; i--)
    { 
      //if (!(*buff) || (*buff == '#')) {number ++; continue;}
      infile1_d.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> eps_1[i] >> prs_1[i] >> ent_1[i] >> nb >> mub;
      delete iss;
    } 

  for(int i=ne_2-1; i>=0; i--)
    { 
      infile2_d.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> eps_2[i] >> prs_2[i] >> ent_2[i] >> nb >> mub;
      delete iss;
    } 

  for(int i=ne_3-1; i>=0; i--)
    { 
      infile3_d.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> eps_3[i] >> prs_3[i] >> ent_3[i] >> nb >> mub;
      delete iss;
    } 

  for(int i=ne_4-1; i>=0; i--)
    { 
      infile4_d.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> eps_4[i] >> prs_4[i] >> ent_4[i] >> nb >> mub;
      delete iss;
    } 

  for(int i=ne_5-1; i>=0; i--)
    { 
      infile5_d.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> eps_5[i] >> prs_5[i] >> ent_5[i] >> nb >> mub;
      delete iss;
    } 

  for(int i=ne_6-1; i>=0; i--)
    { 
      infile6_d.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> eps_6[i] >> prs_6[i] >> ent_6[i] >> nb >> mub;
      delete iss;
    } 

  for(int i=ne_7-1; i>=0; i--)
    { 
      infile7_d.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> eps_7[i] >> prs_7[i] >> ent_7[i] >> nb >> mub;
      delete iss;
    } 

  infile1_t.getline(buff,200);
  infile2_t.getline(buff,200);
  infile3_t.getline(buff,200);
  infile4_t.getline(buff,200);
  infile5_t.getline(buff,200);
  infile6_t.getline(buff,200);
  infile7_t.getline(buff,200);

  infile1_t.getline(buff,200);
  infile2_t.getline(buff,200);
  infile3_t.getline(buff,200);
  infile4_t.getline(buff,200);
  infile5_t.getline(buff,200);
  infile6_t.getline(buff,200);
  infile7_t.getline(buff,200);
  
  for(int i=ne_1-1; i>=0; i--)
    { 
      infile1_t.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> tmp_1[i]  >> mus >> muq;
      delete iss;
    } 

  for(int i=ne_2-1; i>=0; i--)
    { 
      infile2_t.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> tmp_2[i]  >> mus >> muq;
      delete iss;
    } 

  for(int i=ne_3-1; i>=0; i--)
    { 
      infile3_t.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> tmp_3[i]  >> mus >> muq;
      delete iss;
    } 
  for(int i=ne_4-1; i>=0; i--)
    { 
      infile4_t.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> tmp_4[i]  >> mus >> muq;
      delete iss;
    } 
  for(int i=ne_5-1; i>=0; i--)
    { 
      infile5_t.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> tmp_5[i]  >> mus >> muq;
      delete iss;
    } 

  for(int i=ne_6-1; i>=0; i--)
    { 
      infile6_t.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> tmp_6[i]  >> mus >> muq;
      delete iss;
    } 
  for(int i=ne_7-1; i>=0; i--)
    { 
      infile7_t.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> tmp_7[i]  >> mus >> muq;
      delete iss;
    } 


 infile1_d.close();
 infile2_d.close();
 infile3_d.close();
 infile4_d.close();
 infile5_d.close();
 infile6_d.close();
 infile7_d.close();

 infile1_t.close();
 infile2_t.close();
 infile3_t.close();
 infile4_t.close();
 infile5_t.close();
 infile6_t.close();
 infile7_t.close();


   acc = gsl_interp_accel_alloc();

   spline_temp_1 = gsl_spline_alloc (gsl_interp_cspline, ne_1);
   gsl_spline_init (spline_1, eps_1, tmp_1, ne_1);

   spline_temp_2 = gsl_spline_alloc (gsl_interp_cspline, ne_2);
   gsl_spline_init (spline_2, eps_2, tmp_2, ne_2);

   spline_temp_3 = gsl_spline_alloc (gsl_interp_cspline, ne_3);
   gsl_spline_init (spline_3, eps_3, tmp_3, ne_3);

   spline_temp_4 = gsl_spline_alloc (gsl_interp_cspline, ne_4);
   gsl_spline_init (spline_4, eps_4, tmp_4, ne_4);

   spline_temp_5 = gsl_spline_alloc (gsl_interp_cspline, ne_5);
   gsl_spline_init (spline_5, eps_5, tmp_5, ne_5);

   spline_temp_6 = gsl_spline_alloc (gsl_interp_cspline, ne_6);
   gsl_spline_init (spline_6, eps_6, tmp_6, ne_6);

   spline_temp_7 = gsl_spline_alloc (gsl_interp_cspline, ne_7);
   gsl_spline_init (spline_7, eps_7, tmp_7, ne_7);



   spline_prs_1 = gsl_spline_alloc (gsl_interp_cspline, ne_1);
   gsl_spline_init (spline_prs_1, eps_1, prs_1, ne_1);

  spline_prs_2 = gsl_spline_alloc (gsl_interp_cspline, ne_2);
   gsl_spline_init (spline_prs_2, eps_2, prs_2, ne_2);

   spline_prs_3 = gsl_spline_alloc (gsl_interp_cspline, ne_3);
   gsl_spline_init (spline_prs_3, eps_3, prs_3, ne_3);

   spline_prs_4 = gsl_spline_alloc (gsl_interp_cspline, ne_4);
   gsl_spline_init (spline_prs_4, eps_4, prs_4, ne_4);

   spline_prs_5 = gsl_spline_alloc (gsl_interp_cspline, ne_5);
   gsl_spline_init (spline_prs_5, eps_5, prs_5, ne_5);

   spline_prs_6 = gsl_spline_alloc (gsl_interp_cspline, ne_6);
   gsl_spline_init (spline_prs_6, eps_6, prs_6, ne_6);

   spline_prs_7 = gsl_spline_alloc (gsl_interp_cspline, ne_7);
   gsl_spline_init (spline_prs_7, eps_7, prs_7, ne_7);


}


double EoS::pressure( double eg,double _nb, double _nq, double _ns)
{
 if(std::isinf(eg) or std::isnan(eg)) {
   std::cout << "inf energy entered to eos\n"<<std::endl;
  }


        double yi;
        if (eg >= e0_1 && eg < e0_2) { yi = gsl_spline_eval (spline_prs_1, eg, acc);}
        if (eg >= e0_2 && eg < e0_3) { yi = gsl_spline_eval (spline_prs_2, eg, acc);}
        if (eg >= e0_3 && eg < e0_4) { yi = gsl_spline_eval (spline_prs_3, eg, acc);}
        if (eg >= e0_4 && eg < e0_5) { yi = gsl_spline_eval (spline_prs_4, eg, acc);}
        if (eg >= e0_5 && eg < e0_6) { yi = gsl_spline_eval (spline_prs_5, eg, acc);}
        if (eg >= e0_6 && eg < e0_7) { yi = gsl_spline_eval (spline_prs_6, eg, acc);}
        if (eg >= e0_7             ) { yi = gsl_spline_eval (spline_prs_7, eg, acc);}
        if( eg< e0_1 || eg > e0_7  ) { yi = eg/3.0;}

       if(std::isinf(yi) or std::isnan(yi))
         {
          std::cout << "inf pressure returned by eos\n"<<std::endl;
          std::cout<<"energy input is : "<<eg<<std::endl;
        }
        return yi;
}

double EoS::temperature( double eg,double _nb, double _nq, double _ns)
{

        double yi;
        if (eg >= e0_1 && eg < e0_2) { yi = gsl_spline_eval (spline_temp_1, eg, acc);}
        if (eg >= e0_2 && eg < e0_3) { yi = gsl_spline_eval (spline_temp_2, eg, acc);}
        if (eg >= e0_3 && eg < e0_4) { yi = gsl_spline_eval (spline_temp_3, eg, acc);}
        if (eg >= e0_4 && eg < e0_5) { yi = gsl_spline_eval (spline_temp_4, eg, acc);}
        if (eg >= e0_5 && eg < e0_6) { yi = gsl_spline_eval (spline_temp_5, eg, acc);}
        if (eg >= e0_6 && eg < e0_7) { yi = gsl_spline_eval (spline_temp_6, eg, acc);}
        if (eg >= e0_7             ) { yi = gsl_spline_eval (spline_temp_7, eg, acc);}
        return yi;
 
}


double EoS::cs_2( double eg,double _nb, double _nq, double _ns)
{
 if(std::isinf(eg) or std::isnan(eg)) {
   std::cout << "inf energy entered to eos\n"<<std::endl;
  }


        double yi;
        if (eg >= e0_1 && eg < e0_2) { yi = gsl_spline_eval_deriv (spline_prs_1, eg, acc);}
        if (eg >= e0_2 && eg < e0_3) { yi = gsl_spline_eval_deriv (spline_prs_2, eg, acc);}
        if (eg >= e0_3 && eg < e0_4) { yi = gsl_spline_eval_deriv (spline_prs_3, eg, acc);}
        if (eg >= e0_4 && eg < e0_5) { yi = gsl_spline_eval_deriv (spline_prs_4, eg, acc);}
        if (eg >= e0_5 && eg < e0_6) { yi = gsl_spline_eval_deriv (spline_prs_5, eg, acc);}
        if (eg >= e0_6 && eg < e0_7) { yi = gsl_spline_eval_deriv (spline_prs_6, eg, acc);}
        if (eg >= e0_7             ) { yi = gsl_spline_eval_deriv (spline_prs_7, eg, acc);}
        if( eg< e0_1 || eg > e0_7  ) { yi = 1./3.0;}

       if(std::isinf(yi) or std::isnan(yi))
         {
          std::cout << "inf cs2 returned by eos\n"<<std::endl;
        }
        return yi;
}

