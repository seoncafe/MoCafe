  Ò$  F   k820309              17.0        »0Y                                                                                                           
       ../sersic.f90 SERSIC                                                     
       WP DP PI HALFPI                      @                              
       GAMMLN GAMMP INTERP LOCATE                      @                              
       RANDOM_NUMBER                                                      u #RANDOM_KISS32    #RANDOM_KISS64    #RANDOM_KISS_V32    #RANDOM_KISS_V64    #RANDOM_KISS_VV32    #RANDOM_KISS_VV64 	   #RANDOM_KISS_VVV32 
   #RANDOM_KISS_VVV64                      @                              '°                    #NAME    #TAU_FACEON    #RSCALE    #ZSCALE    #RMAX    #ZMAX    #RCENTER                                                                                                ^                                       f              C                                                                                                                                                                                                            
                           ^                    
                                 0.0                                                            
                           ^                    
                       @        6.0                                                            
                           ^                    
                 É?        0.2                                                            
                           ^                    
                       2@        18.0                                                             
                           ^                    
                       2@        18.0                                                   ¨         
                           ^                    
                       @        6.0                      @                              '@             
      #DISKNAME    #RSCALE    #ZSCALE    #RMAX    #ZMAX    #BULGENAME    #SERSIC_INDEX    #REFF    #AXIAL_RATIO    #BULGETODISK                                                                                                                                        r              Cexponential                                                                                                                                                                                                 
                           ^                     
                       @        6.0                                                            
                           ^                     
                     É?        0.2                                                            
                           ^                     
                       2@        18.0                                                            
                           ^                     
                       2@        18.0                                                                                                                                        t              Csersic                                                                                                                                                                                                      
                           ^                     
                       @        4.0                                                   (        
                           ^                     
                    `ff@        2.3                                                   0      	  
                           ^                     
                       à?        0.5                                                   8      
  
                           ^                     
                                 0.0                                                                                                                                                                                                                                                             !     
                
          )       -DTû!	@        3.141592653589793238462643383279502884197                                            "     
                   
                  -DTû!ù?        %         @                               #                    
       #XX $             
                                 $     
      %         @                               %                    
       #A &   #X '             
                                 &     
                
                                 '     
      #         @                                   (                    #NX )   #X *   #Y +   #XNEW ,   #YNEW -             
                                  )                    
                                 *                    
    p          5 O p            5 O p                                   
                                 +                    
    p          5 O p            5 O p                                    
                                 ,     
                                                -     
       #         @                                  .                    #XX /   #N 0   #X 1   #J 2                                            /                    
     p          5 O p            5 O p                                                                     0                                                      1     
                                                  2                                                        3     RANDOM_NUMBER #         @     @                                                #HARVEST 4                                             4     	       #         @     @                                               #HARVEST 5                                             5     
       #         @     @                                                #HARVEST 6                                            6                   	               &                                           #         @     @                                                #HARVEST 7                                            7                   
               &                                           #         @     @                                                #HARVEST 8                                            8                   	               &                   &                                           #         @     @                            	                    #HARVEST 9                                            9                   
               &                   &                                           #         @     @                            
                    #HARVEST :                                            :                   	 	              &                   &                   &                                           #         @     @                                                #HARVEST ;                                            ;                   
               &                   &                   &                                           #         @                                  <                   #SERSIC_CUMULATIVE_3D%BCOEFF =   #M >   #RADIUS ?   #PROFILE @   #RMAX_IN A                                               =                   
                                                    T
W
p          n
         
                  Mî«0':?          n
             
                  m7DS¿]?          n
             
                  û-;Ð¹è?          n
             
                  ÀÉö¾¿          h  p          p          p            p                                                        
                                 >     
                
D@                              ?                   
               &                                                     
D                                @                   
               &                                                     
 @                              A     
      %         @                                 B                    
       #M C   #RMAX D             
  @                              C     
                
 @                              D     
                   fn#fn    ½   P   J  DEFINE      [   J  MATHLIB    h  N   J  RANDOM )   ¶  ê       gen@RANDOM_NUMBER+RANDOM !      £       DUST_TYPE+DEFINE &   C  =  a   DUST_TYPE%NAME+DEFINE ,     §   a   DUST_TYPE%TAU_FACEON+DEFINE (   '  §   a   DUST_TYPE%RSCALE+DEFINE (   Î  §   a   DUST_TYPE%ZSCALE+DEFINE &   u  ¨   a   DUST_TYPE%RMAX+DEFINE &     ¨   a   DUST_TYPE%ZMAX+DEFINE )   Å  §   a   DUST_TYPE%RCENTER+DEFINE #   l  ×       SOURCE_TYPE+DEFINE ,   C	  =  a   SOURCE_TYPE%DISKNAME+DEFINE *   
  §   a   SOURCE_TYPE%RSCALE+DEFINE *   '  §   a   SOURCE_TYPE%ZSCALE+DEFINE (   Î  ¨   a   SOURCE_TYPE%RMAX+DEFINE (   v  ¨   a   SOURCE_TYPE%ZMAX+DEFINE -     =  a   SOURCE_TYPE%BULGENAME+DEFINE 0   [  §   a   SOURCE_TYPE%SERSIC_INDEX+DEFINE (     §   a   SOURCE_TYPE%REFF+DEFINE /   ©  §   a   SOURCE_TYPE%AXIAL_RATIO+DEFINE /   P  §   a   SOURCE_TYPE%BULGETODISK+DEFINE    ÷  p       WP+DEFINE    g  p       DP+DEFINE    ×         PI+DEFINE    p  p       HALFPI+DEFINE    à  X       GAMMLN+MATHLIB "   8  @   a   GAMMLN%XX+MATHLIB    x  ^       GAMMP+MATHLIB     Ö  @   a   GAMMP%A+MATHLIB       @   a   GAMMP%X+MATHLIB    V  r       INTERP+MATHLIB "   È  @   a   INTERP%NX+MATHLIB !     ¤   a   INTERP%X+MATHLIB !   ¬  ¤   a   INTERP%Y+MATHLIB $   P  @   a   INTERP%XNEW+MATHLIB $     @   a   INTERP%YNEW+MATHLIB    Ð  e       LOCATE+MATHLIB "   5  ¤   a   LOCATE%XX+MATHLIB !   Ù  @   a   LOCATE%N+MATHLIB !     @   a   LOCATE%X+MATHLIB !   Y  @   a   LOCATE%J+MATHLIB %     F       RANDOM_NUMBER+RANDOM %   ß  U      RANDOM_KISS32+RANDOM -   4  @   a   RANDOM_KISS32%HARVEST+RANDOM %   t  U      RANDOM_KISS64+RANDOM -   É  @   a   RANDOM_KISS64%HARVEST+RANDOM '   	  U      RANDOM_KISS_V32+RANDOM /   ^     a   RANDOM_KISS_V32%HARVEST+RANDOM '   ê  U      RANDOM_KISS_V64+RANDOM /   ?     a   RANDOM_KISS_V64%HARVEST+RANDOM (   Ë  U      RANDOM_KISS_VV32+RANDOM 0      ¤   a   RANDOM_KISS_VV32%HARVEST+RANDOM (   Ä  U      RANDOM_KISS_VV64+RANDOM 0     ¤   a   RANDOM_KISS_VV64%HARVEST+RANDOM )   ½  U      RANDOM_KISS_VVV32+RANDOM 1     ¼   a   RANDOM_KISS_VVV32%HARVEST+RANDOM )   Î  U      RANDOM_KISS_VVV64+RANDOM 1   #  ¼   a   RANDOM_KISS_VVV64%HARVEST+RANDOM %   ß         SERSIC_CUMULATIVE_3D ,   u   ä     SERSIC_CUMULATIVE_3D%BCOEFF '   Y"  @   a   SERSIC_CUMULATIVE_3D%M ,   "     a   SERSIC_CUMULATIVE_3D%RADIUS -   %#     a   SERSIC_CUMULATIVE_3D%PROFILE -   ±#  @   a   SERSIC_CUMULATIVE_3D%RMAX_IN    ñ#  a       RAND_SERSIC    R$  @   a   RAND_SERSIC%M !   $  @   a   RAND_SERSIC%RMAX 