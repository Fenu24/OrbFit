      program prova 

      ka = 0.D0 

      kd = -(-2.D0*da**3*cos(delta)**2.D0*sin(delta)*ddd+sin(delta)*ddd*&
     &da*dd**2+2.D0*da**2*cos(delta)**2.D0*sin(delta)*dda*dd-sin(delta)*&
     &dda*dd**3-da**5*cos(delta)**3.D0-4.D0*da**3*cos(delta)*dd**2+da**3&
     &*cos(delta)**3.D0*dd**2-2.D0*da*cos(delta)*dd**4)/sqrt(da**2*cos(d&
     &elta)**2.D0+dd**2)**5.D0                                          

      kda = -(2.D0*da**2*cos(delta)**3.D0*ddd-ddd*cos(delta)*dd**2+da**2&
     &*cos(delta)**2.D0*sin(delta)*dd**2-2.D0*sin(delta)*dd**4-3.D0*da*c&
     &os(delta)**3.D0*dda*dd)/sqrt(da**2*cos(delta)**2.D0+dd**2)**5.D0  

      kdd = -(dda*cos(delta)**3.D0*da**2-2.D0*cos(delta)*dda*dd**2-dd*da&
     &**3*sin(delta)*cos(delta)**2.D0+2.D0*da*sin(delta)*dd**3+3.D0*dd*c&
     &os(delta)*ddd*da)/sqrt(da**2*cos(delta)**2.D0+dd**2)**5.D0        

      kdda = -dd*cos(delta)/sqrt(da**2*cos(delta)**2.D0+dd**2)**3.D0 

      kddd = da*cos(delta)/sqrt(da**2*cos(delta)**2.D0+dd**2)**3.D0 


      edota = 0.D0 

      edotd = -da*(da**2*cos(delta)**3.D0*sin(delta)*dda+2.D0*dda*cos(de&
     &lta)*sin(delta)*dd**2-da*dd**3+2.D0*da*dd**3*cos(delta)**2.D0+da**&
     &3*dd*cos(delta)**4.D0-da*cos(delta)*sin(delta)*ddd*dd)/sqrt(da**2*&
     &cos(delta)**2.D0+dd**2)**3.D0                                     

      edotda = -cos(delta)*dd*(-dda*cos(delta)*dd+da**3*cos(delta)**2.D0&
     &*sin(delta)+2.D0*da*dd**2*sin(delta)+da*cos(delta)*ddd)/sqrt(da**2&
     &*cos(delta)**2.D0+dd**2)**3.D0                                    

      edotdd = -da*cos(delta)**2.D0*(-ddd*da+da**3*cos(delta)*sin(delta)&
     &+dda*dd)/sqrt(da**2*cos(delta)**2.D0+dd**2)**3.D0                 

      edotdda = da*cos(delta)**2.D0/sqrt(da**2*cos(delta)**2.D0+dd**2) 

      edotddd = dd/sqrt(da**2*cos(delta)**2.D0+dd**2) 
    END program prova
