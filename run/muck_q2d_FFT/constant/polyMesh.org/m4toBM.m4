// blockMesh :  Block mesh description file
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// ************************************
changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'printf ($1)')])
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])


// atencio Sd manifold trampejat

   meshGenApp blockMesh;
   convertToMeters 1.000;

    define(NB1, 50)
    define(NB2, 10)
    define(NB3, 10)    
    define(NB4, 140)
    define(Ny,  80)
    define(NyB, 36)    
    define(Nyc, 8)
    define(Nz,  100)   

   //******* channel coordinates and mesh ***********
   
          //******* mainfold ***********
   define(Xi, 	0)
   define(Xci,  1)   // 
   define(Xco,  1.2)   // 
   define(Xo,   4) 

   define(Yl, 	2)         // l = left, mirant des del plasma
   define(Ycl,  1.1)
   define(Ycr,  0.9)
   define(Yr, 	0)         // r = right, mirant des del plasma

   define(Zu, 	2)          // u = up
   define(Zd, 	0)          // d = down
   
   vertices
   (
    // INLET B11 **********************************************
    (Xi  Ycr  Zd) vlabel(B11_i_l_d)
    (Xi  Yr  Zd) vlabel(B11_i_r_d)
    (Xi  Ycr  Zu) vlabel(B11_i_l_u)
    (Xi  Yr  Zu) vlabel(B11_i_r_u)
    // OUTLET B11 *********************************************
    (Xci  Ycr  Zd) vlabel(B11_o_l_d)
    (Xci  Yr  Zd) vlabel(B11_o_r_d)
    (Xci  Ycr  Zu) vlabel(B11_o_l_u)
    (Xci  Yr  Zu) vlabel(B11_o_r_u)

    // INLET B12 **********************************************
    (Xi  Ycl  Zd) vlabel(B12_i_l_d)
    (Xi  Ycr  Zd) vlabel(B12_i_r_d)
    (Xi  Ycl  Zu) vlabel(B12_i_l_u)
    (Xi  Ycr  Zu) vlabel(B12_i_r_u)
    // OUTLET B12 *********************************************
    (Xci  Ycl  Zd) vlabel(B12_o_l_d)
    (Xci  Ycr  Zd) vlabel(B12_o_r_d)
    (Xci  Ycl  Zu) vlabel(B12_o_l_u)
    (Xci  Ycr  Zu) vlabel(B12_o_r_u)

    // INLET B13 **********************************************
    (Xi  Yl  Zd) vlabel(B13_i_l_d)
    (Xi  Ycl  Zd) vlabel(B13_i_r_d)
    (Xi  Yl  Zu) vlabel(B13_i_l_u)
    (Xi  Ycl  Zu) vlabel(B13_i_r_u)
    // OUTLET B13 *********************************************
    (Xci  Yl  Zd) vlabel(B13_o_l_d)
    (Xci  Ycl  Zd) vlabel(B13_o_r_d)
    (Xci  Yl  Zu) vlabel(B13_o_l_u)
    (Xci  Ycl  Zu) vlabel(B13_o_r_u)

    // INLET B2 **********************************************
    (Xci  Ycr  Zd) vlabel(B2_i_l_d)
    (Xci  Yr   Zd) vlabel(B2_i_r_d)
    (Xci  Ycr  Zu) vlabel(B2_i_l_u)
    (Xci  Yr   Zu) vlabel(B2_i_r_u)

    // OUTLET B2 **********************************************
    (Xco  Ycr  Zd) vlabel(B2_o_l_d)
    (Xco  Yr   Zd) vlabel(B2_o_r_d)
    (Xco  Ycr  Zu) vlabel(B2_o_l_u)
    (Xco  Yr   Zu) vlabel(B2_o_r_u)

    // INLET B3 **********************************************
    (Xci  Yl    Zd) vlabel(B3_i_l_d)
    (Xci  Ycl   Zd) vlabel(B3_i_r_d)
    (Xci  Yl    Zu) vlabel(B3_i_l_u)
    (Xci  Ycl   Zu) vlabel(B3_i_r_u)

    // OUTLET B3 **********************************************
    (Xco  Yl    Zd) vlabel(B3_o_l_d)
    (Xco  Ycl   Zd) vlabel(B3_o_r_d)
    (Xco  Yl    Zu) vlabel(B3_o_l_u)
    (Xco  Ycl   Zu) vlabel(B3_o_r_u)


    // INLET B41 **********************************************
    (Xco  Ycr  Zd) vlabel(B41_i_l_d)
    (Xco  Yr  Zd) vlabel(B41_i_r_d)
    (Xco  Ycr  Zu) vlabel(B41_i_l_u)
    (Xco  Yr  Zu) vlabel(B41_i_r_u)
    // OUTLET B41 *********************************************
    (Xo  Ycr  Zd) vlabel(B41_o_l_d)
    (Xo  Yr  Zd) vlabel(B41_o_r_d)
    (Xo  Ycr  Zu) vlabel(B41_o_l_u)
    (Xo  Yr  Zu) vlabel(B41_o_r_u)    

    // INLET B42 **********************************************
    (Xco  Ycl  Zd) vlabel(B42_i_l_d)
    (Xco  Ycr  Zd) vlabel(B42_i_r_d)
    (Xco  Ycl  Zu) vlabel(B42_i_l_u)
    (Xco  Ycr  Zu) vlabel(B42_i_r_u)
    // OUTLET B42 *********************************************
    (Xo  Ycl  Zd) vlabel(B42_o_l_d)
    (Xo  Ycr  Zd) vlabel(B42_o_r_d)
    (Xo  Ycl  Zu) vlabel(B42_o_l_u)
    (Xo  Ycr  Zu) vlabel(B42_o_r_u)           
    
    // INLET B43 **********************************************
    (Xco  Yl  Zd) vlabel(B43_i_l_d)
    (Xco  Ycl  Zd) vlabel(B43_i_r_d)
    (Xco  Yl  Zu) vlabel(B43_i_l_u)
    (Xco  Ycl  Zu) vlabel(B43_i_r_u)
    // OUTLET B43 *********************************************
    (Xo  Yl  Zd) vlabel(B43_o_l_d)
    (Xo  Ycl  Zd) vlabel(B43_o_r_d)
    (Xo  Yl  Zu) vlabel(B43_o_l_u)
    (Xo  Ycl  Zu) vlabel(B43_o_r_u)        
       
   );				


   blocks
   (
    //*******************************************
    //*****************************************************************
	// B11
    hex (B11_i_r_d B11_o_r_d B11_o_l_d B11_i_l_d B11_i_r_u B11_o_r_u B11_o_l_u B11_i_l_u)
    (NB1 NyB Nz) simpleGrading ( 1 1 1 )
    // end block B11
    
	// B12
    hex (B12_i_r_d B12_o_r_d B12_o_l_d B12_i_l_d B12_i_r_u B12_o_r_u B12_o_l_u B12_i_l_u)
    (NB1 Nyc Nz) simpleGrading ( 1 1 1 )
    // end block B12    

	// B13
    hex (B13_i_r_d B13_o_r_d B13_o_l_d B13_i_l_d B13_i_r_u B13_o_r_u B13_o_l_u B13_i_l_u)
    (NB1 NyB Nz) simpleGrading ( 1 1 1 )
    // end block B13    
    
    

	// B2
    hex (B2_i_r_d B2_o_r_d B2_o_l_d B2_i_l_d B2_i_r_u B2_o_r_u B2_o_l_u B2_i_l_u)
    (NB2 NyB Nz) simpleGrading ( 1 1 1 )
    // end block B2    

	// B3
    hex (B3_i_r_d B3_o_r_d B3_o_l_d B3_i_l_d B3_i_r_u B3_o_r_u B3_o_l_u B3_i_l_u)
    (NB3 NyB Nz) simpleGrading ( 1 1 1 )
    // end block B3    

	// B41
    hex (B41_i_r_d B41_o_r_d B41_o_l_d B41_i_l_d B41_i_r_u B41_o_r_u B41_o_l_u B41_i_l_u)
    (NB4 NyB Nz) simpleGrading ( 1 1 1 )
    // end block B41
    
	// B42
    hex (B42_i_r_d B42_o_r_d B42_o_l_d B42_i_l_d B42_i_r_u B42_o_r_u B42_o_l_u B42_i_l_u)
    (NB4 Nyc Nz) simpleGrading ( 1 1 1 )
    // end block B42    

	// B43
    hex (B43_i_r_d B43_o_r_d B43_o_l_d B43_i_l_d B43_i_r_u B43_o_r_u B43_o_l_u B43_i_l_u)
    (NB4 NyB Nz) simpleGrading ( 1 1 1 )
    // end block B43  

   );       //end blocks


//   edges
//   (
//     // Xr *************************************
//     arc xr_r1u_r1t xr_r1u_r1b (Xr Yr1arcu Zr1arc0)
//     arc xr_r1u_r1b xr_r1d_r1b (Xr Yr1arc0 Zr1arcb)
//     arc xr_r1d_r1b xr_r1d_r1t (Xr Yr1arcd Zr1arc0)
//     arc xr_r1d_r1t xr_r1u_r1t (Xr Yr1arc0 Zr1arct)
//     );



   patches
   (
      patch xinlet 
      (
	(B11_i_r_d   B11_i_r_u    B11_i_l_u    B11_i_l_d)
    (B12_i_r_d   B12_i_r_u    B12_i_l_u    B12_i_l_d)   	
    (B13_i_r_d   B13_i_r_u    B13_i_l_u    B13_i_l_d)
      )

     patch xoutlet 
      (
	(B41_o_r_d   B41_o_r_u    B41_o_l_u    B41_o_l_d)
	(B42_o_r_d   B42_o_r_u    B42_o_l_u    B42_o_l_d)
	(B43_o_r_d   B43_o_r_u    B43_o_l_u    B43_o_l_d)		
      )

     patch walls
      (
	(B11_i_r_d   B11_o_r_d    B11_o_r_u    B11_i_r_u)
	(B2_i_r_d   B2_o_r_d    B2_o_r_u    B2_i_r_u)	
	(B41_i_r_d   B41_o_r_d    B41_o_r_u    B41_i_r_u)	
	
	(B13_i_l_d   B13_i_l_u    B13_o_l_u    B13_o_l_d)
	(B3_i_l_d   B3_i_l_u    B3_o_l_u    B3_o_l_d)
	(B43_i_l_d   B43_i_l_u    B43_o_l_u    B43_o_l_d)		

    (B11_i_r_u   B11_o_r_u    B11_o_l_u    B11_i_l_u)
    (B12_i_r_u   B12_o_r_u    B12_o_l_u    B12_i_l_u)
    (B13_i_r_u   B13_o_r_u    B13_o_l_u    B13_i_l_u)        
    (B2_i_r_u   B2_o_r_u    B2_o_l_u    B2_i_l_u)
    (B3_i_r_u   B3_o_r_u    B3_o_l_u    B3_i_l_u)
    (B41_i_r_u   B41_o_r_u    B41_o_l_u    B41_i_l_u)   
    (B42_i_r_u   B42_o_r_u    B42_o_l_u    B42_i_l_u)            
    (B43_i_r_u   B43_o_r_u    B43_o_l_u    B43_i_l_u)                     

    (B11_i_r_d   B11_i_l_d    B11_o_l_d    B11_o_r_d)
    (B12_i_r_d   B12_i_l_d    B12_o_l_d    B12_o_r_d)
    (B13_i_r_d   B13_i_l_d    B13_o_l_d    B13_o_r_d)
    (B2_i_r_d   B2_i_l_d    B2_o_l_d    B2_o_r_d)
    (B3_i_r_d   B3_i_l_d    B3_o_l_d    B3_o_r_d)
    (B41_i_r_d   B41_i_l_d    B41_o_l_d    B41_o_r_d)
    (B42_i_r_d   B42_i_l_d    B42_o_l_d    B42_o_r_d)
    (B43_i_r_d   B43_i_l_d    B43_o_l_d    B43_o_r_d)        
        )
    
        patch cleft
        (
    (B3_i_r_d   B3_o_r_d   B3_o_r_u      B3_i_r_u)    
        )    
        patch cright
        (
    (B2_i_l_d   B2_i_l_u    B2_o_l_u    B2_o_l_d)    
        )    
        patch cinlet
        (
    (B12_o_r_u   B12_o_r_d    B12_o_l_d    B12_o_l_u)    
        )
        patch coutlet
        (
    (B42_i_l_d   B42_i_r_d    B42_i_r_u    B42_i_l_u)    
        )    




        patch iil_inter
        (
    (B13_o_r_u   B13_o_r_d    B13_o_l_d    B13_o_l_u)    
        )    
        patch iol_inter
        (
    (B3_i_l_d   B3_i_r_d    B3_i_r_u    B3_i_l_u)    
        )            
        patch oil_inter
        (
    (B3_o_r_u   B3_o_r_d    B3_o_l_d    B3_o_l_u)    
        )    
        patch ool_inter
        (
    (B43_i_l_d   B43_i_r_d    B43_i_r_u    B43_i_l_u)    
        )           
    
        patch iir_inter
        (
    (B11_o_r_u   B11_o_r_d    B11_o_l_d    B11_o_l_u)    
        )    
        patch ior_inter
        (
    (B2_i_l_d   B2_i_r_d    B2_i_r_u    B2_i_l_u)    
        )            
        patch oir_inter
        (
    (B2_o_r_u   B2_o_r_d    B2_o_l_d    B2_o_l_u)    
        )    
        patch oor_inter
        (
    (B41_i_l_d   B41_i_r_d    B41_i_r_u    B41_i_l_u)    
        )       





    
        patch rri_inter
        (
    (B11_i_l_d   B11_i_l_u    B11_o_l_u    B11_o_l_d)    
        )
        patch rli_inter
        (
    (B12_i_r_d   B12_o_r_d   B12_o_r_u      B12_i_r_u)    
        )    
        patch lri_inter
        (
    (B12_i_l_d   B12_i_l_u    B12_o_l_u    B12_o_l_d)    
        )
        patch lli_inter
        (
    (B13_i_r_d   B13_o_r_d   B13_o_r_u      B13_i_r_u)    
        )        
        
        patch rro_inter
        (
    (B41_i_l_d   B41_i_l_u    B41_o_l_u    B41_o_l_d)    
        )
        patch rlo_inter
        (
    (B42_i_r_d   B42_o_r_d   B42_o_r_u      B42_i_r_u)    
        )    
        patch lro_inter
        (
    (B42_i_l_d   B42_i_l_u    B42_o_l_u    B42_o_l_d)    
        )
        patch llo_inter
        (
    (B43_i_r_d   B43_o_r_d   B43_o_r_u      B43_i_r_u)    
        )                
        
        
    
    
     
  );














