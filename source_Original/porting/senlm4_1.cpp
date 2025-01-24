//********1*********1*********1*********1*********1*********1********1**
// Pour ETACHA4
//.....calcula la secc. eficaz total de excitacion desde los
//.....estados ns, np1 y np-1 al estado final nlm, por S.Eikonal
//.....Autor: Cesar A. Ramirez            Fecha: 10-12-2012
//         zt=carga del nucleo blanco
//         zp=carga del nucleo proyectil 
//         ei=energia del estado ligado inicial
//         ef=energia del estado ligado final
//         n,l,m = numeros cuanticos del estado final
//         eta=componente normal del momento transferido
//**********************************************************************
//----- Program to calculte total cross sections of 
//----- monoelectronic atoms excitation, for colision with nude projectiles 
//----- in the SE (Symmetric-Eikonal) aproximation.
//--------------------------------------------------
//......changes in the present version (JPR 02/2012)
//----- zp: Projectile charge (to be excited)
//----- zt: Target charge
//......SC and ASC corrections
//-------------------------------------------------
//----- (n0,l0,m0) quantum numbers of initial state 
//----- (n,l,m) quantum numbers of final state
//----- SUBROUTINES:
//----- * rns, rnp0, rnp1, rnp1m: excitation from ns, np0, np1 np-1 to n'lm
//-----                    respectively with n'lm arbitraries (n'>n)
//----- * intdef: for definite integrals from 0 to infinity
//----- * fac: factorial function
//----- * cg: Clebsch-Gordan coeficients
//----- * armonic: spherical harmonic function
//----- * plgndr: Legendre polinomials
//----- * gamac: gama function
//----- * cf21d, chyper, hypgfx: hypergeometric functions F21(ca,cb,cc,cx)
//**********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SEnlm(zp8,zt8,E8) {
      // implicit double*8[a-b][d-h][o-z],complex*16[c];
      double*4 SeSE,StSE;
      char kk*2;
      }
       if(dabs[csum] > preci*dabs[cint]){
      goto L39;
      }
 else {
      } //continue;
      } ;
      return;
      }
//***********************************************************************