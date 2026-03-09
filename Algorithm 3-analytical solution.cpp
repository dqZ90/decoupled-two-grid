// two-level grad-div stabilization method
// The natural convection equation based on Ossen linearization
// Decoupled algorithm 3

//relations between coarse and fine mesh
// mr=1: H=h^2/3, h=n^3;
// mr=2: H=h^1/2, h=n^2;
// mr=3: H=m*h, h=n*G;
// mr=4: H=h^1/2, h=2^n;
int[int] mra(4);
mra[0]=1;mra[1]=2;mra[2]=3;mra[3]=4;
int mri=1, mrf=1; // mri and mrf are the initial and final arry mra respectively;
int step=2; // meshing refining step
int N0=6;   // meshing initial size
int N1=6;     // corresponding to mr=1
int N2=14;     // corresponding to mr=2
int N3=12;     // corresponding to mr=3
int N4=7;     // corresponding to mr=4
int G=10;    // mesh triangulation parameter: 1/h=G*n;
int m=2;   // H=m*h; H/m=h^2

/// nonlinear iteration scheme
//NOS=1: Newton scheme
//NOS=2: Ossen  scheme
//NOS=3: Stokes scheme 
int NOSi=2, NOSf=2; //NOSi and NOSf are the initial and final values of NOS, respectively;
int K=5000; //admisible maximum number of spatial nonlinear iterations 

// for different Ra
real[int] Raa(9);
Raa[0]=1;Raa[1]=10;Raa[2]=100;Raa[3]=1000;Raa[4]=10000; Raa[5]=100000;Raa[6]=1000000;Raa[7]=10000000;Raa[8]=10000000;
int Rai=3, Raf=3;//nui and nuf are the initial and final arry Paa respectively;
int Rait;
real Pr=1,  kap=1, Ra; // parameter in the equation
real Sigma=100;//grad-div stabilization parameter


//stopping criterion for nonlinear iteration:
//  ||U^(n+1)-U^n||_0/||U^(n+1)||_0}<TOL;


//typies of mesh on nonoverlapping subregions
//case=1: structrual mesh;
//case=2: unstructrual mesh;
int cas=1;

int NOS,NOSit,mr,mrit,nuit,Ni,Nf;
real epsi=1.0e-10; //parameter in the stablized term for computing
real TOL=1.0e-6;    // the error tolerance of the iteration 

real UL2norm,UH1norm,UH2norm, TL2norm,TH1norm,TH2norm,PL2norm, PH1norm; // norms of the exact solution on whole domain
real eherr,ehnorm; 
real errUL2,errUH1,errPL2, errTL2,errTH1, errwhole,errsum,olderrUL2,olderrUH1,olderrPL2, olderrTL2,olderrTH1;
real olderrwhole,olderrsum,UL2rate,UH1rate,PL2rate, TL2rate,TH1rate, wholerate,sumrate;
real iterrUL2, iterrUH1, lnorm, iterr;
int k,hk,Hk; 
real CPU,CPU0;

real uc=0.1;
real pc=10;

// exact steady solutions
func u1tmpx=x^2-2*x^3+x^4;//=x^2*(x-1)^2
func u1tmpy=y-3*y^2+2*y^3;//=y*(y-1)*(2*y-1)
func u2tmpx=x-3*x^2+2*x^3; //=x*(x-1)*(2*x-1)
func u2tmpy=y^2-2*y^3+y^4; //=y^2*(y-1)^2

func u1exact=uc*u1tmpx*u1tmpy;
func u2exact=-uc*u2tmpy*u2tmpx;
func Texact=u1exact+u2exact;
func pexact=pc*(6*x^5+6*y^5-2);

func u1exactdx=uc*(2*x-6*x^2+4*x^3)*u1tmpy;
func u1exactdxx=uc*(2-12*x+12*x^2)*u1tmpy;
func u1exactdy=uc*u1tmpx*(1-6*y+6*y^2);
func u1exactdyy=uc*u1tmpx*(-6+12*y);
func u1exactdxdy=uc*(2*x-6*x^2+4*x^3)*(1-6*y+6*y^2);

func u2exactdx=-uc*u2tmpy*(1-6*x+6*x^2);
func u2exactdxx=-uc*u2tmpy*(-6+12*x);
func u2exactdy=-uc*(2*y-6*y^2+4*y^3)*u2tmpx;
func u2exactdyy=-uc*(2-12*y+12*y^2)*u2tmpx;
func u2exactdxdy=-uc*(2*y-6*y^2+4*y^3)*(1-6*x+6*x^2);

func Texactdx=u1exactdx+u2exactdx;
func Texactdy=u1exactdy+u2exactdy;
func Texactdxx=u1exactdxx+u2exactdxx;
func Texactdyy=u1exactdyy+u2exactdyy;
func Texactdxdy=u1exactdxdy+u2exactdxdy;

func pexactdx=pc*6*5*x^4;
func pexactdy=pc*6*5*y^4;

//f=-Pr*delta*U+(U\cdot\nabla)U+gradp-Pr*Ra*e_j*T
func f1=-Pr*(u1exactdxx+u1exactdyy)+u1exact*u1exactdx+u2exact*u1exactdy+pexactdx;
func f2=-Pr*(u2exactdxx+u2exactdyy)+u1exact*u2exactdx+u2exact*u2exactdy+pexactdy-Pr*Ra*Texact;

//g=-kap*delta*T+(u\cdot\naba)T
func g=-kap*(Texactdxx+Texactdyy)+u1exact*Texactdx+u2exact*Texactdy;
//compute the norms of the exact solutions 
{  
   mesh Th0=square(100,100);
   UL2norm=sqrt( int2d(Th0)((u1exact^2+u2exact^2)));
   UH1norm=sqrt( int2d(Th0)((u1exactdx^2+u1exactdy^2+u2exactdx^2+u2exactdy^2)));
   UH2norm=sqrt( int2d(Th0)((u1exactdxx^2+2*u1exactdxdy^2+u1exactdyy^2
                           +u2exactdxx^2+2*u2exactdxdy^2+u2exactdyy^2)));
   TL2norm=sqrt( int2d(Th0)(Texact^2));
   TH1norm=sqrt( int2d(Th0)((Texactdx^2+Texactdy^2)));
   TH2norm=sqrt( int2d(Th0)((Texactdxx^2+2*Texactdxdy^2+Texactdyy^2)));
   PL2norm=sqrt( int2d(Th0)(pexact^2));
   PH1norm=sqrt( int2d(Th0)(pexactdx^2+pexactdy^2));
}


for(int nosj=NOSi;nosj<=NOSf;nosj++){
   NOS=nosj;  
   for(mrit=mri;mrit<=mrf;mrit++){ 
      mr=mra[mrit]; 
      if(mr==1) {Ni=N0;Nf=N1;}
      if(mr==2) {Ni=N0;Nf=N2;}
      if(mr==3) {Ni=N0;Nf=N3;}
      if(mr==4) {Ni=N0;Nf=N4;}
   
          for(Rait=Rai;Rait<=Raf;Rait++){
             Ra=Raa[Rait];

             string E ="P2-P1-P2",S="NC-2level-polynomial";
             string R;
             if(mr==1) R=",h=n^3,H^3=h^2,";
             if(mr==2) R=",h=n^2,H^2=h,";
             if(mr==3) R=",h=n*"+G+",H="+m+"h,";
             if(mr==4) R=",h=2^n,H^2=2h,";
              
           
             ofstream ffcout("decoupled-Pr="+Pr+"-Ra="+Ra+"-kap="+kap+"-"+E+"(Sigma="+Sigma+").txt");

   ffcout<<S+".edp results\n";
   ffcout<<"Decoupled two-grid method for NC equations\n";   
   ffcout<<"Example:  polynomial function exact solution\n";
   ffcout<<"Omega=[0,1]*[0,1]"<<endl; 
   ffcout<<"Sigma="<<Sigma<<"   ,Ra="<<Ra<<endl;
   ffcout<<"b(u,v,w)=0.5*((u,grad)v,w)-0.5*((u,grad)w,v)\n";
  
  
   if(mr==3)ffcout<<"H= "<<m<<"h"<<endl;   
   ffcout<<E+" element, the parameter in 'epsi*p*q' for stabilized computing epsi= "<<epsi<<endl; 
  
   ffcout<<"The nonlinear iteration method  is Oseen ";
 
   
   ffcout<<"\nStopping criterion for nolinear iteration :\n";   
   ffcout<<"||U^(n+1)-U^n||_0,Omega/||U^(n+1)||_0,Omega}<"<<TOL<<endl<<endl;           
  

for (int n=Ni;n<=Nf;n=n+step){
     if(mr==1){hk=n*n*n; Hk=n*n;} //H=h^2/3, h=n^-3       
     if(mr==2){hk=n*n;  Hk=2*n;} //H=h^1/2, h=n^-2 
     if(mr==3){hk=n*G; Hk=hk/m;} //H=m*h,h=1/n*G         
     real h=1.0/hk, H=1.0/Hk;
   
     ffcout<<"\n1/h="<<hk<<",    1/H="<<Hk<<endl;

     //The first step: solving the problems on a coarse mesh
    
     real CPU0=clock();
     //coarse grid nonlinear problem  
     mesh TH=square(Hk,Hk);   
     fespace Uh(TH,P2); Uh u1,u2,v1,v2, oldu1,oldu2; // the velocity space
     fespace Ph(TH,P1); Ph p,q; // the pressure space 
     fespace Teh(TH,P2); Teh T,Phi,oldT; // the pressure space 

     //Stokes problem for the initial value of nonlinear iteration    
     problem Stokes([u1,u2,p, T],[v1,v2,q, Phi]) =
            int2d(TH)(Pr*(dx(u1)*dx(v1)+dy(u1)*dy(v1)+ dx(u2)*dx(v2)+ dy(u2)*dy(v2))//a(U,V)                     
            -p*(dx(v1)+dy(v2))  //d(V,p)
            +q*(dx(u1)+dy(u2))    //d(U,Q)
            -Pr*Ra*T*v2
            +epsi*p*q           // stabilized term for computing  
            +Sigma*(dx(u1)*dx(v1)+dx(u1)*dy(v2)+dy(u2)*dx(v1)+dy(u2)*dy(v2))
            +kap*(dx(T)*dx(Phi)+dy(T)*dy(Phi)))
            -int2d(TH)(f1*v1+f2*v2+g*Phi)   //right hand side
            +on(1,2,3,4,u1=0,u2=0, T=0);   //boundary condition     

     //scheme2: b(u,v,w)=0.5*((u,grad)v,w)-0.5*((u,grad)w,v);
     problem NC01([u1,u2,p, T],[v1,v2,q,Phi]) =
            int2d(TH)(Pr*(dx(u1)*dx(v1)+dy(u1)*dy(v1)+dx(u2)*dx(v2)+dy(u2)*dy(v2))//a(U,V)
            +0.5*(oldu1*dx(u1)*v1+oldu2*dy(u1)*v1+oldu1*dx(u2)*v2+oldu2*dy(u2)*v2) //Trilinear term b(U,U,V)
            -0.5*(oldu1*dx(v1)*u1+oldu2*dy(v1)*u1+oldu1*dx(v2)*u2+oldu2*dy(v2)*u2)
			+0.5*(u1*dx(oldu1)*v1+u2*dy(oldu1)*v1+u1*dx(oldu2)*v2+u2*dy(oldu2)*v2)
			-0.5*(u1*dx(v1)*oldu1+u2*dy(v1)*oldu1+u1*dx(v2)*oldu2+u2*dy(v2)*oldu2) //Trilinear term b(U,oldU,V)
            -p*(dx(v1)+dy(v2))  //d(V,p)
            +q*(dx(u1)+dy(u2))    //d(U,Q)
            -Pr*Ra*T*v2
            +epsi*p*q             // stabilized term for computing 
            +Sigma*(dx(u1)*dx(v1)+dx(u1)*dy(v2)+dy(u2)*dx(v1)+dy(u2)*dy(v2))  //stabilization trem
            +kap*(dx(T)*dx(Phi)+dy(T)*dy(Phi))
            +0.5*(oldu1*dx(T)*Phi+oldu2*dy(T)*Phi) //Trilinear term b(u,T,phi)
            -0.5*(oldu1*dx(Phi)*T+oldu2*dy(Phi)*T)
			+0.5*(u1*dx(oldT)*Phi+u2*dy(oldT)*Phi)
			-0.5*(u1*dx(Phi)*oldT+u2*dy(Phi)*oldT))
            -int2d(TH)(f1*v1+f2*v2+g*Phi
			+0.5*(oldu1*dx(oldu1)*v1+oldu2*dy(oldu1)*v1+oldu1*dx(oldu2)*v2+oldu2*dy(oldu2)*v2)
			-0.5*(oldu1*dx(v1)*oldu1+oldu2*dy(v1)*oldu1+oldu1*dx(v2)*oldu2+oldu2*dy(v2)*oldu2) //Trilinear term b(oldu,oldu,v)
			+0.5*(oldu1*dx(oldT)*Phi+oldu2*dy(oldT)*Phi)
			-0.5*(oldu1*dx(Phi)*oldT+oldu2*dy(Phi)*oldT))   //right hand side
            +on(1,2,3,4,u1=0,u2=0, T=0);   //boundary condition

   
    int it=0;
    real ite=1;        
    Stokes;    //solving Stokes problem to initiate the iteration 
    while(ite>TOL && it<K){
        it=it+1;
        oldu1=u1;
        oldu2=u2;
        oldT=T;
        NC01;
           
    //compute the successious error                
        iterrUL2 = sqrt( int2d(TH)((u1-oldu1)^2+(u2-oldu2)^2 ));
        lnorm=sqrt( int2d(TH)((u1)^2+(u2)^2 ));
        if(lnorm<1.0e-20) lnorm=1.0e-20;
        ite=iterrUL2/lnorm;
        }
     
    if(it==K && ite>TOL){
        ffcout<<"\n      Nonlinear iterating "<<K<<" didn't satisfy the stopping condition on coarse grid !\n";
        break;
        }
        else {ffcout<<"\n      Number of  nonlinear iterations on coarse grid= "<<it<<endl;}

    //The second step: fine grid nonlinear problem
          
    mesh Th1=square(hk,hk); 
    fespace Uh01(Th1,P2); Uh01 eh11,eh12,oldeh11,oldeh12,vh11,vh12,u11=u1,u12=u2; // the velocity space
    fespace Ph01(Th1,P1); Ph01 rh1,qh1,p1=p;     // the pressure space
    fespace Teh01(Th1,P2); Teh01 eT1, phi1,TT1=T;     // the pressure space

    problem NCOssen1([eh11,eh12,rh1],[vh11,vh12,qh1]) =
            int2d(Th1)(Pr*(dx(eh11)*dx(vh11)+dy(eh11)*dy(vh11)+dx(eh12)*dx(vh12)+dy(eh12)*dy(vh12))
			            +0.5*(u11*dx(eh11)*vh11+u12*dy(eh11)*vh11+u11*dx(eh12)*vh12+u12*dy(eh12)*vh12) //Trilinear term b(U,U,V)
                        -0.5*(u11*dx(vh11)*eh11+u12*dy(vh11)*eh11+u11*dx(vh12)*eh12+u12*dy(vh12)*eh12)
						-rh1*(dx(vh11)+dy(vh12))
						+qh1*(dx(eh11)+dy(eh12))
						+epsi*rh1*qh1
						+Sigma*(dx(eh11)*dx(vh11)+dx(eh11)*dy(vh12)+dy(eh12)*dx(vh11)+dy(eh12)*dy(vh12)))
			-int2d(Th1)(f1*vh11+f2*vh12+Pr*Ra*TT1*vh12)
			+on(1,2,3,4,eh11=0,eh12=0);
    NCOssen1;
    oldeh11=eh11;
    oldeh12=eh12;
    problem NCOssen2([eT1],[phi1]) =
		  int2d(Th1)(kap*(dx(eT1)*dx(phi1)+dy(eT1)*dy(phi1))
		             +0.5*(oldeh11*dx(eT1)*phi1+oldeh12*dy(eT1)*phi1) //Trilinear term b(u,T,phi)
                     -0.5*(oldeh11*dx(phi1)*eT1+oldeh12*dy(phi1)*eT1))
		  -int2d(Th1)(g*phi1)
		  +on(1,2,3,4,eT1=0);

    NCOssen2;	
    //compute the errors 
	CPU=clock()-CPU0;     
    errUH1=sqrt(int2d(Th1)((dx(eh11)-u1exactdx)^2
                      +(dy(eh11)-u1exactdy)^2
                      +(dx(eh12)-u2exactdx)^2
                      +(dy(eh12)-u2exactdy)^2));
    errPL2=sqrt(int2d(Th1)((pexact - rh1)^2));     
    errTH1=sqrt(int2d(Th1)((dx(eT1)-Texactdx)^2
                      +(dy(eT1)-Texactdy)^2));         
 
    real CPU1=clock()-CPU0;
    ffcout<<"\n      CPU time (grid+computing)="<<CPU<<"\n " ;
    ffcout<<"\n      CPU time (grid+computing+error computing)= "<<CPU1<<"\n " ;
    ffcout<<"      errUH1 = "<<errUH1<<endl;
    ffcout<<"      errPL2 = "<<errPL2<<endl;
    ffcout<<"      errTH1 = "<<errTH1<<endl;
 
      
if(n>Ni){ 
    real n1,n2;
    if(mr==1) n1=(n-step)*(n-step)*(n-step); //H^3=h^2,h=n^-3;             
    if(mr==2) n1=(n-step)*(n-step);      //H=h^2,h=n^-2; 
    if(mr==3) n1=G*(n-step);           //H=mh,h=G*n;
    if(mr==4) n1=2^(n-step);           //H=h^1/2, h=2^-n;     
    n2=hk;
    real hh=log(n2)-log(n1);          
    UH1rate=(log(olderrUH1)-log(errUH1))/hh;
    PL2rate=(log(olderrPL2)-log(errPL2))/hh;
    TH1rate=(log(olderrTH1)-log(errTH1))/hh;
                
          
    ffcout<<"      U_H1 convergence rate= "<<UH1rate<<endl;
    ffcout<<"      P_L2 convergence rate= "<<PL2rate<<endl;
    ffcout<<"      T_H1 convergence rate= "<<TH1rate<<endl;	
	
}
    olderrUH1=errUH1;
    olderrPL2=errPL2;
    olderrTH1=errTH1;


}
    ffcout<<"U_H1 norm="<<UH1norm<<endl;
    ffcout<<"P_L2 norm="<<PL2norm<<endl;
    ffcout<<"T_H1 norm="<<TH1norm<<endl;    
    
}  
}
}
 



