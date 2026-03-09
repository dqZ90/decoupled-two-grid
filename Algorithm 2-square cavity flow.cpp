//A two-level grad-div stabilizations method based on the P2-P1-P2 element
//Global coarse grid and local fine grid linearied correection method
//example:Buoyancy-driven square cavity flow
//Algorithm 2

//relations between coarse and fine mesh
// mr=1: H=h^2/3, h=n^3;
// mr=2: H=h^1/2, h=n^2;
// mr=3: H=m*h, h=n*G;
// mr=4: H=h^1/2, h=2^n;
int[int] mra(4);
mra[0]=1;mra[1]=2;mra[2]=3;mra[3]=4;
int mri=2, mrf=2; // mri and mrf are the initial and final arry mra respectively;
int step=1; // meshing refining step
int N0=1;   // meshing initial size
int N1=1;     // corresponding to mr=1
int N2=1;     // corresponding to mr=2
int N3=1;     // corresponding to mr=3
int N4=1;     // corresponding to mr=3
int G=162;    // mesh triangulation parameter: 1/h=G*n;
real m=2;   // H=m*h;

// nonlinear iteration method on coarse grid
//NOS=1: Newton method
//NOS=2: Ossen method
//NOS=3: Simple method

int[int] NOSa(3);
NOSa[0]=1; NOSa[1]=2; NOSa[2]=3;    
int NOSi=0, NOSf=0; //NOSi and NOSf are the initial and final arry NOSa respectively; 
int M=5000;  //Maximum nonlinear ieteration steps

//stopping criterion for nonlinear iteration on coarse mesh:
// sign=0: max{||U^(n+1)-U^n||_0,(Omega)/||U^(n+1)||_0,(Omega)}<TOL;
// sign=1: |logh|^{1/2}||grad(U^(n+1)-U^n)||_0,Omega||U^(n+1)-U^n||_0,Omega}<"+c+"H^2";
// sign=2: ||U^(n+1)-U^n||_{0,Omega}<"+c+"H^2";

int sign=0;

// linearized methods on fine grid
//LNOS=1: Newton-linearized method
//LNOS=2: Ossen-linearized  method
//LNOS=3: Simple-linearized  method
int[int] LNOSa(3);
LNOSa[0]=1; LNOSa[1]=2; LNOSa[2]=3;
int LNOSi=1, LNOSf=1; //LNOSi and LNOSf are the initial and final arry LNOSa respectively; 
//NOte: only Newton-linearized method is available

// scheme for the trilinear term
//scheme=1: b(u,v,w)=((u,grad)v,w)
//scheme=2: b(u,v,w)=0.5*((u,grad)v,w)-0.5*((u,grad)w,v);
//scheme=3: b(u,v,w)=((u,grad)v,w)+0.5*((div u v),w)
int scheme=2;

// for different Ra
real[int] Raa(4);
Raa[0]=1000;Raa[1]=10000;Raa[2]=100000;Raa[3]=1000000;
int Rai=0, Raf=0;//Rai and Raf are the initial and final arry Raa respectively;
int Rait;
real Pr=0.71,  kap=1, Ra; // parameters in the equation
real Sigma=10;  //stabilization parameter

int NOS,NOSit,LNOS,LNOSit,mr,mrit,nuit,Ni,Nf;
real TOL=1.0e-6;   //stopping criterion for nolinear iteration
real epsi=1.0e-8; //parameter in the stablized term for computing
//int drw=30; //Once a drw iteration steps to plot

real UL2norm,UH1norm,UH2norm,TL2norm,Thnorm,TH2norm, PL2norm, PH1norm; // norms of the exact solution on whole domain
real errUL2,errUH1,errTL2,errTh,errPL2,errwhole,errsum,olderrUL2,olderrUH1,olderrTL2,olderrTh;
real olderrPL2,olderrwhole,olderrsum,UL2rate,UH1rate,TL2rate,Thrate,PL2rate,wholerate,sumrate;
int iit,iite,iiterUL2,iiterrUH1,iiterr,localiiterr,eherr,ehnorm; 
real iterrUL2,lnorm;
real n1,n2;
int k,hk,Hk,xk,yk;
real x0,x1,yo,y1;
real CPU,CPU0;

for(NOSit=NOSi;NOSit<=NOSf; NOSit++){
   NOS=NOSa[NOSit];
   for(LNOSit=LNOSi;LNOSit<=LNOSf; LNOSit++){
     LNOS=LNOSa[LNOSit];
      for(mrit=mri;mrit<=mrf;mrit++){ 
         mr=mra[mrit]; 
      if(mr==1) {Ni=N0;Nf=N1;}
      if(mr==2) {Ni=N0;Nf=N2;}
      if(mr==3) {Ni=N0;Nf=N3;}
      if(mr==4) {Ni=N0;Nf=N4;}   
          for(Rait=Rai;Rait<=Raf;Rait++){
             Ra=Raa[Rait];

             string E ="P2-P1-P2",S="NC-Two-grid-Cavity-Coupled";
             string R,I,L;
             if(mr==1) R=",h=n^3,H^3=h^2,";
             if(mr==2) R=",h=n^2,H^2=h,";
             if(mr==3) R=",n*"+G+",H="+m+"h,"; 
			 if(mr==4) R=",h=2^n,H^2=2h,";
             if(NOS==1) I="-Newton-Method";
             if(NOS==2) I="-Ossen-Method";
             if(NOS==3) I="-Simple-Method";
             if(LNOS==1) L="-LNewton";
             if(LNOS==2) L="-LOssen";
             if(LNOS==3) L="-LSimple";             
                           
   ofstream ffcout("Coupled-Cavity-Ra="+Ra+"(Sigma="+Sigma+").text");

   ffcout<<S+".edp results\n";
   ffcout<<"Global coarse grid and local fine grid linearied correection method based on two-grid mesh for NC equations\n";   
   ffcout<<"Example:  Buoyancy-driven square cavity flow\n";
   ffcout<<"Omega=[0,1]*[0,1]"<<endl;
   if(mr==1)ffcout<<"H= "<<m<<"h"<<endl;
   if(mr==2)ffcout<<"H= "<<m<<"h "<<endl;   
   if(mr==3)ffcout<<"H= "<<m<<"h"<<endl;
   if(mr==4)ffcout<<"H= "<<m<<"h"<<endl;   
   ffcout<<E+" element, the parameter in 'epsi*p*q' for stabilized computing epsi= "<<epsi<<endl; 
   ffcout<<" Ra="<<Ra<<endl;
   ffcout<<" Sigma="<<Sigma<<endl;
   ffcout<<"The nonlinear iteration method on coarse method is: ";
   if(NOS==1)ffcout<<"Newton iteration\n";
   if(NOS==2)ffcout<<"Ossen iteration\n";   
   if(NOS==3)ffcout<<"Simple iteration\n";
   ffcout<<"The local linearized method on fine method is: ";
   if(LNOS==1)ffcout<<"Newton method\n";
   if(LNOS==2)ffcout<<"Ossen method\n";   
   if(LNOS==3)ffcout<<"Simple method\n";

   if(scheme==1) ffcout<<"  b(u,v,w)=((u,grad)v,w)\n";
   if(scheme==2) ffcout<<"  b(u,v,w)=0.5*((u,grad)v,w)-0.5*((u,grad)w,v)\n";
   if(scheme==3) ffcout<<"  b(u,v,w)=((u,grad)v,w)+0.5*((div u v),w)\n";

   if(sign==0)ffcout<<"          max{||U^(n+1)-U^n||_0/||U^(n+1)||_0}<"<<TOL<<endl<<endl;
   if(sign==1)ffcout<<"          max{||U^(n+1)-U^n||_1/||U^(n+1)||_1}<"<<TOL<<endl<<endl;
   
       
     for (int n=Ni;n<=Nf;n=n+step){
     if(mr==1){k=n*n*n;Hk=n*n;} //H=h^2/3, h=n^-3       
     if(mr==2){k=n*n;Hk=n;} //H=h^1/2, h=n^-2 
     if(mr==3){k=n*G;Hk=k/m;} //H=m*h,h=1/n*G
     if(mr==4){k=2^n;Hk=sqrt(k);} //H=h^1/2,h=2^-n         
     
	 real h=1.0/k, H=1.0/Hk;
     ffcout<<"\n 1/h="<<k<<",1/H="<<Hk<<endl;
	 
     border a(t=0,1){x=t; y=0; label=1;};
     border b(t=0,1){x=1; y=t; label=2;};
     border c(t=1,0){x=t; y=1; label=3;};
     border d(t=1,0){x=0; y=t; label=4;};

		
     real CPU0=clock();
     //coarse grid nolinear problem
     mesh TH=square(Hk,Hk);
     fespace Uh01(TH,P2); Uh01 u1,u2,oldu1,oldu2,v1,v2; // the velocity space
     fespace Ph01(TH,P1); Ph01 p,q;     // the pressure space
     fespace Teh01(TH,P2); Teh01 T,Phi,oldT; // the temperature space
 
     //Stokes problem for the initial value of nonlinear iteration of NC 
 
     problem NCStokes([u1,u2,p,T],[v1,v2,q,Phi])=
        int2d(TH)(Pr*(dx(u1)*dx(v1)+dy(u1)*dy(v1)+dx(u2)*dx(v2)+dy(u2)*dy(v2))   //Pr*a(U,V)
	             -p*(dx(v1)+dy(v2))  //d(v,p)
			     +q*(dx(u1)+dy(u2))  //d(u,q)
			     -Pr*Ra*T*v2        //Pr*Ra*(e_jT,v)
			     +epsi*p*q          //stabilized term for computing
			     +Sigma*(dx(u1)*dx(v1)+dx(u1)*dy(v2)+dy(u2)*dx(v1)+dy(u2)*dy(v2))  //grad-div stabilization term
	             +kap*(dx(T)*dx(Phi)+dy(T)*dy(Phi)))  //Kap*(T,Phi)
		         +on(1,3,u1=0,u2=0)  //boundary condition
			     +on(2,u1=0,u2=0,T=0)
			     +on(4,u1=0,u2=0,T=1);

     //NCOssen method
     problem NCNewton([u1,u2,p,T],[v1,v2,q,Phi])=
        int2d(TH)(Pr*(dx(u1)*dx(v1)+dy(u1)*dy(v1)+dx(u2)*dx(v2)+dy(u2)*dy(v2))  //Pr*a(U,V)
                 +0.5*(oldu1*dx(u1)*v1+oldu2*dy(u1)*v1+oldu1*dx(u2)*v2+oldu2*dy(u2)*v2) //Trilinear term b(oldU,U,V)
                 -0.5*(oldu1*dx(v1)*u1+oldu2*dy(v1)*u1+oldu1*dx(v2)*u2+oldu2*dy(v2)*u2)
				 +0.5*(u1*dx(oldu1)*v1+u2*dy(oldu1)*v1+u1*dx(oldu2)*v2+u2*dy(oldu2)*v2)
				 -0.5*(u1*dx(v1)*oldu1+u2*dy(v1)*oldu1+u1*dx(v2)*oldu2+u2*dy(v2)*oldu2) //Trilinear term b(U,oldU,V)
				 -p*(dx(v1)+dy(v2))  //d(v,p)
			     +q*(dx(u1)+dy(u2))  //d(u,q)
			     -Pr*Ra*T*v2        //Pr*Ra*(e_jT,v)
			     +epsi*p*q          //stabilized term for computing
				 +Sigma*(dx(u1)*dx(v1)+dx(u1)*dy(v2)+dy(u2)*dx(v1)+dy(u2)*dy(v2))  //grad-div stabilization
				 +kap*(dx(T)*dx(Phi)+dy(T)*dy(Phi))  //Kap*(T,Phi)
				 +0.5*(oldu1*dx(T)*Phi+oldu2*dy(T)*Phi) //Trilinear term b(oldu,T,phi)
                 -0.5*(oldu1*dx(Phi)*T+oldu2*dy(Phi)*T)
				 +0.5*(u1*dx(oldT)*Phi+u2*dy(oldT)*Phi)
				 -0.5*(u1*dx(Phi)*oldT+u2*dy(Phi)*oldT)) //Trilinear term b(u,oldT,phi)
       -int2d(TH)(
	             +0.5*(oldu1*dx(oldu1)*v1+oldu2*dy(oldu1)*v1+oldu1*dx(oldu2)*v2+oldu2*dy(oldu2)*v2)
			     -0.5*(oldu1*dx(v1)*oldu1+oldu2*dy(v1)*oldu1+oldu1*dx(v2)*oldu2+oldu2*dy(v2)*oldu2) //Trilinear term b(oldu,oldu,v)
			     +0.5*(oldu1*dx(oldT)*Phi+oldu2*dy(oldT)*Phi)
			     -0.5*(oldu1*dx(Phi)*oldT+oldu2*dy(Phi)*oldT))		//Trilinear term b(oldu,oldT,phi)				
	   +on(1,3,u1=0,u2=0)  //boundary condition
	   +on(2,u1=0,u2=0,T=0)
	   +on(4,u1=0,u2=0,T=1);

       //solve the NC equations
       int it=0;
       real ite=1.0;
       NCStokes; //solving Stokes problem to initiate the iteration
       while(ite>TOL && it<M){
            it=it+1;
            oldu1=u1;
            oldu2=u2;
			oldT=T;
            NCNewton;	
        
      //compute the successious error
      iterrUL2 = sqrt( int2d(TH)((u1-oldu1)^2+(u2-oldu2)^2 ));
      lnorm=sqrt( int2d(TH)((u1)^2+(u2)^2 ));
      //avoid divided by zero
      if(lnorm<1.0e-10) lnorm=1.0e-10;
             ite=iterrUL2/lnorm;
        }  
      if(it>=M && ite>TOL){
      ffcout<<"  Nonlinear iterating "<<M<<" didn't satisfy the stopping condition on subdomain1 !\n";
        }
	  else{ffcout<<"\n Number of nonlinear iterations on coarse grid="<<it<<endl;}
	  ffcout<<"  CPU time for coarse grid problem="<<CPU0<<endl;

      //fine grid nonlinear problem
      mesh Th=square(k,k);
      fespace Uh1(Th,P2);Uh1 eh11,eh12,vh11,vh12,u11=u1,u12=u2;//the velocity space
      fespace Ph1(Th,P1);Ph1 rh1,qh1;//the pressure space
      fespace The1(Th,P2);The1 eT1,phi1,TT1=T;//the pressure space

      problem Newton([eh11,eh12,rh1,eT1],[vh11,vh12,qh1,phi1])=
           int2d(Th)(Pr*(dx(eh11)*dx(vh11)+dy(eh11)*dy(vh11)+dx(eh12)*dx(vh12)+dy(eh12)*dy(vh12)) //a(U,V)
           +0.5*(u11*dx(eh11)*vh11+u12*dy(eh11)*vh11+u11*dx(eh12)*vh12+u12*dy(eh12)*vh12) //Trilinear term b(U,U,V)
           -0.5*(u11*dx(vh11)*eh11+u12*dy(vh11)*eh11+u11*dx(vh12)*eh12+u12*dy(vh12)*eh12)
           -rh1*(dx(vh11)+dy(vh12))  //d(V,p)
           +qh1*(dx(eh11)+dy(eh12))    //d(U,Q)
           -Pr*Ra*eT1*vh12
           +epsi*rh1*qh1             // stabilized term for computing
		   +Sigma*(dx(eh11)*dx(vh11)+dx(eh11)*dy(vh12)+dy(eh12)*dx(vh11)+dy(eh12)*dy(vh12))  //stabilization term
           +kap*(dx(eT1)*dx(phi1)+dy(eT1)*dy(phi1))
           +0.5*(u11*dx(eT1)*phi1+u12*dy(eT1)*phi1) //Trilinear term b(u,T,phi)
           -0.5*(u11*dx(phi1)*eT1+u12*dy(phi1)*eT1))
           +on(1,3,eh11=0,eh12=0)//boundary condition
           +on(2,eh11=0,eh12=0,eT1=0)
           +on(4,eh11=0,eh12=0,eT1=1);

      Newton;
      CPU=clock()-CPU0;
      real CPU1=clock()-CPU0;

     real intNua=int1d(Th,4)(-dx(eT1));
     ffcout<<" Nusselt average value on x=1: "<<intNua<<endl;
     ffcout<<"\n 1/h="<<k<< ", 1/H="<<Hk<<endl;
     ffcout<<"\n       CPU time (grid+ computing )= " <<CPU<< endl;
     ffcout<<"        CPU time (grid+ computing+ error comuting)= "<<CPU1<<endl;
 if(M>0){
     ffcout<<"\n      number nonlinear iterations on whole domain= "<<it<<endl; 
}   

     //for tecplot
     int i,j;
     ofstream fu1("U1-Ra="+Ra+"-"+I+"-"+E+".dat");
     ofstream fp1("P1-Ra="+Ra+"-"+I+"-"+E+".dat");
     ofstream ft1("T1-Ra="+Ra+"-"+I+"-"+E+".dat");

         fu1<<"Variables="<<"X"<<","<<"Y"<<","<<"u"<<","<<"v"<<endl;
         fu1<<"Zone"<<"   "<<"N="<<Th.nv<<","<<"E="<<Th.nt<<","<<"F=FEPOINT,ET=TRIANGLE"<<endl;
         
		 fp1<<"Variables="<<"X"<<","<<"Y"<<","<<"p"<<endl;
         fp1<<"Zone"<<"   "<<"N="<<Th.nv<<","<<"E="<<Th.nt<<","<<"F=FEPOINT,ET=TRIANGLE"<<endl;
         
		 ft1<<"Variables="<<"X"<<","<<"Y"<<","<<"T"<<endl;
         ft1<<"Zone"<<"   "<<"N="<<Th.nv<<","<<"E="<<Th.nt<<","<<"F=FEPOINT,ET=TRIANGLE"<<endl;
                
         for(i=0;i<Th.nv;i++){
             fu1<<Th(i).x<<"   "<<Th(i).y<<"   "<<eh11(Th(i).x,Th(i).y)<<"   "<<eh12(Th(i).x,Th(i).y)<<endl;
             fp1<<Th(i).x<<"   "<<Th(i).y<<"   "<<rh1(Th(i).x,Th(i).y)<<endl;
             ft1<<Th(i).x<<"   "<<Th(i).y<<"   "<<eT1(Th(i).x,Th(i).y)<<endl;
         } 
         for(i=0;i<Th.nt;i++){
            for(j=0;j<3;j++){
               fu1<<Th[i][j]+1<<"  ";
               fp1<<Th[i][j]+1<<"  ";
               ft1<<Th[i][j]+1<<"  ";
            }
            if(j==3){fu1<<endl; fp1<<endl; ft1<<endl;}
			}
       //For the numerical results of velocity along vertical and horizontal lines through geometric center of the cavity
	    ofstream fx("fine grid-X= "+Ra+"-"+I+"-"+E+".text");
        ofstream fuu1("fine grid-U1= "+Ra+"-"+I+"-"+E+".text");
        ofstream fuu2("fine grid-U2= "+Ra+"-"+I+"-"+E+".text");
		ofstream ftty("fine grid-Ty= "+Ra+"-"+I+"-"+E+".text");
		ofstream fNuh0("fine grid-Nuh0= "+Ra+"-"+I+"-"+E+".text");
		ofstream fNuh1("fine grid-Nuh1= "+Ra+"-"+I+"-"+E+".text");
		fx<<"The coordinations of vertical and horizontal lines through geometric center of the cavity:\n\n";
        fuu1<<"u1 velocity along vertical line through geometric center of the cavity:\n\n";
        fuu2<<"u2 velocity along horizontal line through geometric center of the cavity:\n\n";
        ftty<<"T along horizontal line through geometric center of the cavity:\n\n";
		fNuh0<<"Nu=-\partil T/\partial x  along x=0 (hot wall)  of the cavity:\n\n";
		fNuh1<<"Nu=-\partil T/\partial x  along x=1 (cool wall)  of the cavity:\n\n";
		real u1max=-100, u2max=-100,Numax=-100, Numin=1000;
		real u1t,u2t,Nuth0,Nutc,Nuth1;
		for(i=0;i<=k;i++){
			real xy=i*h;
			fx<<xy <<" ";
			     u1t=eh11(0.5,xy);
                 u2t=eh12(xy,0.5);
				 Nuth0=-dx(eT1)(0,xy);
				 Nuth1=-dx(eT1)(1,xy);				 
			fuu1<<u1t<<" ";
	        fuu2<<u2t<<" ";
			fNuh0<<Nuth0<<" ";
			fNuh1<<Nuth1<<" ";
			ftty<<eT1(xy,0.5)<<" ";
			if (u1t>u1max) u1max=u1t;
            if (u2t>u2max) u2max=u2t;
			if (Nuth0>Numax) Numax=Nuth0;
			if (Nuth0<Numin) Numin=Nuth0; 
		}
	 ffcout<<"      Max u1 along  vertical centre line:  "<<u1max<<endl;
     ffcout<<"      Max u2 along  horizontal centre line:  "<<u2max<<endl;
     ffcout<<"      Max -\partil T/\partial x on x=0:  "<<Numax<<endl;
     ffcout<<"      Min -\partil T/\partial x on x=0:  "<<Numin<<endl;
	 }
		  }
	  }
   }
}