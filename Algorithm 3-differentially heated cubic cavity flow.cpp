//A two-level grad-div stabilizations method based on the P2-P1-P2 element
//Global coarse grid and local fine grid linearied correection method
//example: Differentially heated cubic cavity flow
//Algorithm 3

//relations between coarse and fine mesh
// mr=1: H=h^2/3, h=n^3;
// mr=2: H=h^1/2, h=n^2;
// mr=3: H=m*h, h=n*G;
// mr=4: H=h^1/2, h=2^n;
int[int] mra(4);
mra[0]=1;mra[1]=2;mra[2]=3;mra[3]=4;
int mri=2, mrf=2; // mri and mrf are the initial and final arry mra respectively;
int step=2; // meshing refining step
int N0=1;   // meshing initial size
int N1=1;     // corresponding to mr=1
int N2=1;     // corresponding to mr=2
int N3=1;     // corresponding to mr=3
int N4=1;     // corresponding to mr=3
int G=12;    // mesh triangulation parameter: 1/h=G*n;
real m=2;   // H=m*h;

//Support 3D mesh generation and 3D visualization functions
load "msh3" 
load "medit"  // dynamics load tools for 3d.

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
real Pr=1,  kap=1, Ra; // parameters in the equation
real Sigma=1;  //stabilization parameter

int NOS,NOSit,LNOS,LNOSit,mr,mrit,nuit,Ni,Nf;
real TOL=1.0e-8;   //stopping criterion for nolinear iteration
real epsi=1.0e-10; //parameter in the stablized term for computing 
real iterrUL2,lnorm;
//real n1,n2;
int k,hk,Hk,xk,yk;

func mesh3 Cube(int[int] & NN,real[int,int] &BB ,int[int,int] & L){
    
// first build the 6 faces of the hex.
    real x0=BB(0,0),x1=BB(0,1);
    real y0=BB(1,0),y1=BB(1,1);
    real z0=BB(2,0),z1=BB(2,1);
    int nx=NN[0],ny=NN[1],nz=NN[2];
    mesh Thx = square(nx,ny,[x0+(x1-x0)*x,y0+(y1-y0)*y]);
    int[int] rup=[0,L(2,1)], rdown=[0,L(2,0)],
     rmid=[1,L(1,0), 2,L(0,1), 3, L(1,1), 4, L(0,0) ];
     mesh3 Th=buildlayers(Thx,nz, zbound=[z0,z1],labelmid=rmid,
                           labelup = rup,labeldown = rdown);
     return Th;
    }

int [int,int] LL=[[1,2],[3,4],[5,6]];


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

            string E ="P23d-P1-P23d",S="NC-Two-grid-decoupled-cavity";
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
                           
   ofstream ffcout("Decoupled-Cavity-Ra="+Ra+"(Sigma="+Sigma+").text");

   ffcout<<S+".edp results\n";
   ffcout<<"Global coarse grid and local fine grid linearied correection method based on two-grid mesh for NC equations\n";   
   ffcout<<"Example: Differentially heated cubic cavity flow\n";
   ffcout<<"Omega=[0,1]*[0,1]*[0,1]"<<endl;
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
	              
  //coarse grid nolinear problem
  real [int,int] BB=[[0,1],[0,1],[0,1]]; // bounding box of the solution domain
  int[int] NN=[Hk,Hk,Hk]; // the Ramber of step in each direction
  mesh3 TH=Cube(NN,BB,LL); 
  fespace Uh(TH,P23d); Uh u1,u2,u3,v1,v2,v3,oldu1,oldu2,oldu3; // the velocity space
  fespace Ph(TH,P13d); Ph p,q; // the pressure space
  fespace Teh(TH,P23d);Teh T,Phi,oldT;//the temperature space

  //Stokes problem for the initial value of nonlinear iteration of NC 
  problem NCStokes([u1,u2,u3,p,T],[v1,v2,v3,q,Phi]) =
            int3d(TH)(Pr*(dx(u1)*dx(v1)+dy(u1)*dy(v1)+dz(u1)*dz(v1)
			             +dx(u2)*dx(v2)+dy(u2)*dy(v2)+dz(u2)*dz(v2)
                         +dx(u3)*dx(v3)+dy(u3)*dy(v3)+dz(u3)*dz(v3))//a(U,V)
                      -p*(dx(v1)+dy(v2)+dz(v3))  //d(V,p)
                      +q*(dx(u1)+dy(u2)+dz(u3))    //d(U,Q)
					  -Pr*Ra*T*v3
                      +epsi*p*q           // stabilized term for computing
					  +Sigma*(dx(u1)*dx(v1)+dx(u1)*dy(v2)+dx(u1)*dz(v3)
						     +dy(u2)*dx(v1)+dy(u2)*dy(v2)+dy(u2)*dz(v3)
						     +dz(u3)*dx(v1)+dz(u3)*dy(v2)+dz(u3)*dz(v3))//stabilization term
					  +kap*(dx(T)*dx(Phi)+dy(T)*dy(Phi)+dz(T)*dz(Phi)))
			+on(1,u1=0,u2=0,u3=0,T=-1.0/2.0)
			+on(2,u1=0,u2=0,u3=0,T=1.0/2.0)
			+on(3,4,5,6,u1=0,u2=0,u3=0);   //boundary condition  

  //NCNewton method
  problem NCNewton([u1,u2,u3,p,T],[v1,v2,v3,q,Phi]) =
            int3d(TH)(Pr*(dx(u1)*dx(v1)+dy(u1)*dy(v1)+dz(u1)*dz(v1)
			             +dx(u2)*dx(v2)+dy(u2)*dy(v2)+dz(u2)*dz(v2)
                         +dx(u3)*dx(v3)+dy(u3)*dy(v3)+dz(u3)*dz(v3))
					  +0.5*(oldu1*dx(u1)*v1+oldu2*dy(u1)*v1+oldu3*dz(u1)*v1
					       +oldu1*dx(u2)*v2+oldu2*dy(u2)*v2+oldu3*dz(u2)*v2
						   +oldu1*dx(u3)*v3+oldu2*dy(u3)*v3+oldu3*dz(u3)*v3)
                      -0.5*(oldu1*dx(v1)*u1+oldu2*dy(v1)*u1+oldu3*dz(v1)*u1
					       +oldu1*dx(v2)*u2+oldu2*dy(v2)*u2+oldu3*dz(v2)*u2
						   +oldu1*dx(v3)*u3+oldu2*dy(v3)*u3+oldu3*dz(v3)*u3)    //Trilinear term a(oldU,U,V)
                      +0.5*(u1*dx(oldu1)*v1+u2*dy(oldu1)*v1+u3*dz(oldu1)*v1
					       +u1*dx(oldu2)*v2+u2*dy(oldu2)*v2+u3*dz(oldu2)*v2
						   +u1*dx(oldu3)*v3+u2*dy(oldu3)*v3+u3*dz(oldu3)*v3)
					  -0.5*(u1*dx(v1)*oldu1+u2*dy(v1)*oldu1+u3*dz(v1)*oldu1
					       +u1*dx(v2)*oldu2+u2*dy(v2)*oldu2+u3*dz(v2)*oldu2
						   +u1*dx(v3)*oldu3+u2*dy(v3)*oldu3+u3*dz(v3)*oldu3)    //Trilinear term a(U,oldU,V)
					  +Sigma*(dx(u1)*dx(v1)+dx(u1)*dy(v2)+dx(u1)*dz(v3)
						     +dy(u2)*dx(v1)+dy(u2)*dy(v2)+dy(u2)*dz(v3)
							 +dz(u3)*dx(v1)+dz(u3)*dy(v2)+dz(u3)*dz(v3))//stabilization term
					  -p*(dx(v1)+dy(v2)+dz(v3))
					  +q*(dx(u1)+dy(u2)+dz(u3))
					  -Pr*Ra*T*v3
					  +epsi*p*q
					  +kap*(dx(T)*dx(Phi)+dy(T)*dy(Phi)+dz(T)*dz(Phi))
					  +0.5*(oldu1*dx(T)*Phi+oldu2*dy(T)*Phi+oldu3*dz(T)*Phi)
					  -0.5*(oldu1*dx(Phi)*T+oldu2*dy(Phi)*T+oldu3*dz(Phi)*T)   //Trilinear term a(oldU,T,Phi)
					  +0.5*(u1*dx(oldT)*Phi+u2*dy(oldT)*Phi+u3*dz(oldT)*Phi)
					  -0.5*(u1*dx(Phi)*oldT+u2*dy(Phi)*oldT+u3*dz(Phi)*oldT))   //Trilinear term a(U,oldT,Phi)
		   -int3d(TH)(
		              +0.5*(oldu1*dx(oldu1)*v1+oldu2*dy(oldu1)*v1+oldu3*dz(oldu1)*v1
					       +oldu1*dx(oldu2)*v2+oldu2*dy(oldu2)*v2+oldu3*dz(oldu2)*v2
						   +oldu1*dx(oldu3)*v3+oldu2*dy(oldu3)*v3+oldu3*dz(oldu3)*v3)
					  -0.5*(oldu1*dx(v1)*oldu1+oldu2*dy(v1)*oldu1+oldu3*dz(v1)*oldu1
					       +oldu1*dx(v2)*oldu2+oldu2*dy(v2)*oldu2+oldu3*dz(v2)*oldu2
						   +oldu1*dx(v3)*oldu3+oldu2*dy(v3)*oldu3+oldu3*dz(v3)*oldu3)     //Trilinear term a(oldU,oldU,V)
					  +0.5*(oldu1*dx(oldT)*Phi+oldu2*dy(oldT)*Phi+oldu3*dz(oldT)*Phi)
					  -0.5*(oldu1*dx(Phi)*oldT+oldu2*dy(Phi)*oldT+oldu3*dz(Phi)*oldT))  //Trilinear term a(oldU,oldT,Phi)
		   +on(1,u1=0,u2=0,u3=0,T=-1.0/2.0)
		   +on(2,u1=0,u2=0,u3=0,T=1.0/2.0)
		   +on(3,4,5,6,u1=0,u2=0,u3=0);   //boundary condition 

  //solve the NC equations
  int it=0;
  real ite=1.0;

  NCStokes; //solving Stokes problem to initiate the iteration
while(ite>TOL && it<M){
        it=it+1;
        oldu1=u1;
        oldu2=u2;
		oldu3=u3;
		oldT=T;
   NCNewton;	
        
   //compute the successious error
   iterrUL2 = sqrt( int3d(TH)((u1-oldu1)^2+(u2-oldu2)^2 +(u3-oldu3)^2));
   lnorm=sqrt( int3d(TH)((u1)^2+(u2)^2+(u3)^2));
   //avoid divided by zero
   if(lnorm<1.0e-20) lnorm=1.0e-20;
        ite=iterrUL2/lnorm;
        }
   if(it>=M && ite>TOL){
       ffcout<<"  Nonlinear iterating "<<M<<" didn't satisfy the stopping condition on subdomain1 !\n";
        }
		else{ffcout<<"\n Number of nonlinear iterations on coarse grid="<<it<<endl;}	

  //fine grid nonlinear problem
  real [int,int] B1=[[0,1],[0,1],[0,1]]; // bounding box of the solution domain
  int[int] NN1=[k,k,k]; // the Ramber of step in each direction    
  mesh3 Th01=Cube(NN1,B1,LL);
  fespace Uh1(Th01,P23d); Uh1 eh11,eh12,eh13,oldeh11,oldeh12,oldeh13,vh11,vh12,vh13,u11=u1,u12=u2,u13=u3; // the velocity space
  fespace Ph1(Th01,P13d); Ph1 rh1,qh1,p1=p;  //mh1; // the pressure space 
  fespace Teh1(Th01,P23d);Teh1 eT1,phi1,TT1=T; //the temperature space
  problem NCOssen01([eh11,eh12,eh13,rh1],[vh11,vh12,vh13,qh1])=
    int3d(Th01)(Pr*(dx(eh11)*dx(vh11)+dy(eh11)*dy(vh11)+dz(eh11)*dz(vh11)
                   +dx(eh12)*dx(vh12)+dy(eh12)*dy(vh12)+dz(eh12)*dz(vh12)
		           +dx(eh13)*dx(vh13)+dy(eh13)*dy(vh13)+dz(eh13)*dz(vh13)) //a(U,V)
              +0.5*(u11*dx(eh11)*vh11+u12*dy(eh11)*vh11+u13*dz(eh11)*vh11
		           +u11*dx(eh12)*vh12+u12*dy(eh12)*vh12+u13*dz(eh12)*vh12
				   +u11*dx(eh13)*vh13+u12*dy(eh13)*vh13+u13*dz(eh13)*vh13) //Trilinear term b(U,U,V)
              -0.5*(u11*dx(vh11)*eh11+u12*dy(vh11)*eh11+u13*dz(vh11)*eh11
                   +u11*dx(vh12)*eh12+u12*dy(vh12)*eh12+u13*dz(vh12)*eh12
				   +u11*dx(vh13)*eh13+u12*dy(vh13)*eh13+u13*dz(vh13)*eh13)
              -rh1*(dx(vh11)+dy(vh12)+dz(vh13))  //d(V,p)
              +qh1*(dx(eh11)+dy(eh12)+dz(eh13))    //d(U,Q)
              +epsi*rh1*qh1             // stabilized term for computing
		      +Sigma*(dx(eh11)*dx(vh11)+dx(eh11)*dy(vh12)+dx(eh11)*dz(vh13)
		             +dy(eh12)*dx(vh11)+dy(eh12)*dy(vh12)+dy(eh12)*dz(vh13)
				     +dz(eh13)*dx(vh11)+dz(eh13)*dy(vh12)+dz(eh13)*dz(vh13)))  //stabilization term
   -int3d(Th01)(Pr*Ra*TT1*vh13)
   +on(1,2,3,4,5,6,eh11=0,eh12=0,eh13=0);
  NCOssen01;
  oldeh11=eh11;
  oldeh12=eh12;
  oldeh13=eh13;
  problem NCOssen02([eT1],[phi1])=
     int3d(Th01)(kap*(dx(eT1)*dx(phi1)+dy(eT1)*dy(phi1)+dz(eT1)*dz(phi1))
	           +0.5*(oldeh11*dx(eT1)*phi1+oldeh12*dy(eT1)*phi1+oldeh13*dz(eT1)*phi1) //Trilinear term b(u,T,phi)
               -0.5*(oldeh11*dx(phi1)*eT1+oldeh12*dy(phi1)*eT1+oldeh13*dz(phi1)*eT1))
    +on(1,eT1=-1.0/2.0)
	+on(2,eT1=1.0/2.0);
  NCOssen02;

  real intNua=int2d(Th01,2)(dx(eT1));
  ffcout<<" Nusselt average value on x=1: "<<intNua<<endl;
  real eh11max=eh11[].max;
  real eh12max=eh12[].max;
  real eh13max=eh13[].max;
  ffcout<<" eh11 max= "<<eh11max<<endl;
  ffcout<<" eh12 max= "<<eh12max<<endl;
  ffcout<<" eh13 max= "<<eh13max<<endl;

  int gg=22;//22
  real [int,int] B2=[[0,1],[0,1],[0,1]]; // bounding box of the solution domain
  int[int] NN2=[gg,gg,gg]; // the Ramber of step in each direction    
  mesh3 Th02=Cube(NN2,B2,LL);
  fespace Uh2(Th02,P23d); Uh2 eh21=eh11,eh22=eh12,eh23=eh13; // the velocity space
  fespace Ph2(Th02,P13d); Ph2 rh2=rh1;  //mh1; // the pressure space 
  fespace Teh2(Th02,P23d);Teh2 eT2=eT1; //the temperature space

  //for tecplot
  int i,j;
  ofstream fu("U-Ra="+Ra+"-"+I+"-"+E+".dat");
  ofstream fp("P-Ra="+Ra+"-"+I+"-"+E+".dat");
  ofstream ft("T-Ra="+Ra+"-"+I+"-"+E+".dat");

         fu<<"Variables="<<"X"<<","<<"Y"<<","<<"Z"<<","<<"u"<<","<<"v"<<","<<"w"<<endl;
         fu<<"Zone"<<"   "<<"N="<<Th02.nv<<","<<"E="<<Th02.nt<<","<<"F=FEPOINT,ET=TETRAHEDRON"<<endl;
         
		 fp<<"Variables="<<"X"<<","<<"Y"<<","<<"Z"<<","<<"p"<<endl;
         fp<<"Zone"<<"   "<<"N="<<Th02.nv<<","<<"E="<<Th02.nt<<","<<"F=FEPOINT,ET=TETRAHEDRON"<<endl;
         
		 ft<<"Variables="<<"X"<<","<<"Y"<<","<<"Z"<<","<<"T"<<endl;
         ft<<"Zone"<<"   "<<"N="<<Th02.nv<<","<<"E="<<Th02.nt<<","<<"F=FEPOINT,ET=TETRAHEDRON"<<endl;
                
         for(i=0;i<Th02.nv;i++){
             fu<<Th02(i).x<<"   "<<Th02(i).y<<"   "<<Th02(i).z<<"    "<<eh21(Th02(i).x,Th02(i).y,Th02(i).z)<<"   "<<eh22(Th02(i).x,Th02(i).y,Th02(i).z)<<"   "<<eh23(Th02(i).x,Th02(i).y,Th02(i).z)<<endl;
             fp<<Th02(i).x<<"   "<<Th02(i).y<<"   "<<Th02(i).z<<"    "<<rh2(Th02(i).x,Th02(i).y,Th02(i).z)<<endl;
             ft<<Th02(i).x<<"   "<<Th02(i).y<<"   "<<Th02(i).z<<"    "<<eT2(Th02(i).x,Th02(i).y,Th02(i).z)<<endl;
         } 
         for(i=0;i<Th02.nt;i++){
            for(j=0;j<4;j++){
               fu<<Th02[i][j]+1<<"  ";
               fp<<Th02[i][j]+1<<"  ";
               ft<<Th02[i][j]+1<<"  ";
            }
			if(j=4){
				fu<<endl; 
			    fp<<endl; 
			    ft<<endl;
			}
			} 
	 }
		  }
	  }
   }
}