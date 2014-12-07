      program main
      implicit Double Precision (A-H,O-Z)  
      PARAMETER (PI=3.141592653589)
      
      DIMENSION Bx1(1:3), Bx2(1:3), Bx3(1:3), Bx4(1:3)
      
      
      parameter (NII=40000)
      dimension IR(NII),JR(NII),KR(NII),LR(NII),MR(NII),NR(NII)      
      dimension X1(NII),Y1(NII),Z1(NII),X2(NII),Y2(NII),Z2(NII)
      dimension X3(NII),Y3(NII),Z3(NII),X4(NII),Y4(NII),Z4(NII)
     
      parameter (NIJ=13) 
      
      OPEN(8,FILE='setup.inp',STATUS='old')
      READ(8,*) s
      READ(8,*) NN
      READ(8,*) r0  
      READ(8,*) xi  
      READ(8,*) num
      READ(8,*) ssc
      CLOSE(8)
      print*, s, NN, r0, xi, ssc
**********************************************************************

         WRITE(*,*) 'Reading space points from position.dat'
         I=1
        OPEN(20,FILE='position.dat',STATUS='old')
  992     CONTINUE

         READ(20,*,END=1002)  IR(I),JR(I),KR(I),LR(I),NNN,
     &           X1(I),Y1(I),Z1(I),X2(I),Y2(I),Z2(I), 
     &           X3(I),Y3(I),Z3(I),X4(I),Y4(I),Z4(I)             
C     &         AX1(I,J,K,L,M,N),Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4 


         I=I+1
         GOTO 992
 1002    CONTINUE
         Close(20)
         
         Ndata=I-1
      print*, 'datalength', Ndata  
      
      If (Ndata .gt. NII) then
        print*, 'alert, Ndata>NII'
        stop
      endif  
         
      DO 140 I=1,Ndata        
     
     

      sum1 =0.
      sum2 =0.
          
          
      Bx1(1) = X1(I)*ssc
      Bx1(2) = Y1(I)*ssc
      Bx1(3) = Z1(I)*ssc
      
      Bx2(1) = X2(I)*ssc
      Bx2(2) = Y2(I)*ssc
      Bx2(3) = Z2(I)*ssc 
      
      Bx3(1) = X3(I)*ssc
      Bx3(2) = Y3(I)*ssc
      Bx3(3) = Z3(I)*ssc 
      
      Bx4(1) = X4(I)*ssc
      Bx4(2) = Y4(I)*ssc
      Bx4(3) = Z4(I)*ssc
          
         call Fun4(Bx1,Bx2,Bx3,Bx4,s,xi,r0,NN,sum1,sum2) 
         
         Print*, I, sum1,sum2
         
          write(99,*) IR(I),JR(I),KR(I),LR(I), sum1, sum2 
     
 140  CONTINUE		        
 
		  
      end program
      
***************************************************************************

*************************************************************************
**  THIS FUNCTION GIVES THE 4 POINT SPACE CORRELATION OF SIGMA FIELD   **
*************************************************************************
      SUBROUTINE Fun4(Bx1,Bx2,Bx3,Bx4,s,xi,r0,NN,sum1,sum2) !s is the spacing of dz, xi is the correlation length, r0 is the cutoff of integration.
      Implicit Double Precision (A-H,O-Z)
      PARAMETER (PI=3.141592653589)
      DIMENSION Bx1(1:3), Bx2(1:3), Bx3(1:3), Bx4(1:3)      
      
        
          sum1 = 0.
          sum2 = 0.
          
        
          do i=-NN, NN
           do j=-NN, NN
            do k=-NN, NN    

           r1=sqrt((s*i)**2+(s*j)**2+(s*k)**2)		   
           r2=sqrt((s*i-(Bx2(1)-Bx1(1)))**2+(s*j-(Bx2(2)-Bx1(2)))**2
     &         +(s*k-(Bx2(3)-Bx1(3)))**2)	
           r3=sqrt((s*i-(Bx3(1)-Bx1(1)))**2+(s*j-(Bx3(2)-Bx1(2)))**2
     &        +(s*k-(Bx3(3)-Bx1(3)))**2)	
           r4=sqrt((s*i-(Bx4(1)-Bx1(1)))**2+(s*j-(Bx4(2)-Bx1(2)))**2
     &        +(s*k-(Bx4(3)-Bx1(3)))**2)		
		   		   
           if(r1 .lt. 0.001)then
              r1=r0
           end if
           if(r2 .lt. 0.001)then
              r2=r0
           end if
           if(r3 .lt. 0.001)then
              r3=r0
           end if
           if(r4 .lt. 0.001)then
              r4=r0
           end if
		   
           sum1=sum1+exp(-(r1+r2+r3+r4)/xi)/(r1*r2*r3*r4)*s**3		   
		   
            end do
           end do
          end do  
***************************************************************************
!$OMP parallel do reduction(+:sum2) private(r1,r2,r3,r4,r5)
       do i2=-NN, NN 
        do j2=-NN, NN
         do k2=-NN, NN     

          do i=-NN, NN
           do j=-NN, NN
            do k=-NN, NN    

           r1=sqrt((s*i)**2+(s*j)**2+(s*k)**2)	
           r2=sqrt((s*i-(Bx2(1)-Bx1(1)))**2+(s*j-(Bx2(2)-Bx1(2)))**2
     &         +(s*k-(Bx2(3)-Bx1(3)))**2)	
           r3=sqrt((s*i2-(Bx3(1)-Bx1(1)))**2+(s*j2-(Bx3(2)-Bx1(2)))**2
     &        +(s*k2-(Bx3(3)-Bx1(3)))**2)	
           r4=sqrt((s*i2-(Bx4(1)-Bx1(1)))**2+(s*j2-(Bx4(2)-Bx1(2)))**2
     &        +(s*k2-(Bx4(3)-Bx1(3)))**2)	           	
           r5=sqrt((s*i-s*i2)**2+(s*j-s*j2)**2+(s*k-s*k2)**2)	
		   		   
           if(r1 .lt. 0.001)then
              r1=r0
           end if
           if(r2 .lt. 0.001)then
              r2=r0
           end if
           if(r3 .lt. 0.001)then
              r3=r0
           end if
           if(r4 .lt. 0.001)then
              r4=r0
           end if
           if(r5 .lt. 0.001)then
              r5=r0
           end if
		   
           sum2=sum2+exp(-(r1+r2+r3+r4+r5)/xi)/(r1*r2*r3*r4*r5)*s**6		   
		   
            end do
           end do
          end do        

         end do
        end do	
c       print*, i2, sum2		   
       end do
!$OMP END parallel do

          
          sum1 = sum1/(4*PI)**4 
          sum2 = sum2/(4*PI)**5 
      
c          print*, sum1, sum2
      
       Return
       End      
      

************************************************************************
**  THIS FUNCTION GIVES THE 3 POINT SPACE CORRELATION OF SIGMA FIELD  **
************************************************************************
      Function Fun3(x1,x2,x3,s,xi,r0)      ! s is the spacing of dz, xi is the correlation length, r0 is the cutoff of integration.
      Implicit Double Precision (A-H, O-Z)
      PARAMETER (PI=3.141592653589)
      DIMENSION x1(1:3), x2(1:3), x3(1:3)   
      
!        print*, s,xi,r0
      
         NN =20
         sum1 = 0.0
         
         do i=-NN, NN
           do j=-NN, NN
            do k=-NN, NN    

           r1=sqrt((s*i)**2+(s*j)**2+(s*k)**2)	
           r2=sqrt((s*i-(x2(1)-x1(1)))**2+(s*j-(x2(2)-x1(2)))**2
     &         +(s*k-(x2(3)-x1(3)))**2)	
           r3=sqrt((s*i-(x3(1)-x1(1)))**2+(s*j-(x3(2)-x1(2)))**2
     &        +(s*k-(x3(3)-x1(3)))**2)	
		   		   
           if(r1 .lt. 0.001)then
              r1=r0
           end if
           if(r2 .lt. 0.001)then
              r2=r0
           end if
           if(r3 .lt. 0.001)then
              r3=r0
           end if
		   
           sum1=sum1+exp(-(r1+r2+r3)/xi)/(r1*r2*r3)*s**3
          		   
            end do
           end do
          end do        
          Fun3 = sum1/(4*PI)**3 
!          print*, Fun3
      return    
      end      
***************************************************************************
***************************************************************************      
      
