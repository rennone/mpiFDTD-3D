#include "ntff3D.h"
#include <math.h>
#include "field.h"

// 係数と最低限必要なインデックスをまとめて求める為のマクロ
#define COEFF_AND_INDICES(i, j, k)                      \
  double r2x = i+0.5-cx, r2y = j+0.5-cy, r2z = k+0.5-cz; \
  double dot = r1x*r2x + r1y*r2y + r1z*r2z;               \
  dcomplex coef = cexp( I*k_s*dot);                       \
  int w = field_index(i, j, k);

//必ず COEFF_AND_INDICESの後に記述する
//前後の面に必要なE,Hを求める
#define EH_IN_XY(w, ex, ey, hx, hy)                                      \
  int w_rt = field_right(w);                                            \
  int w_tp = field_top(w);                                              \
  int w_ft = field_front(w);                                            \
  dcomplex ex = -0.5*( Ex[w] + Ex[w_rt] );                              \
  dcomplex ey =  0.5*( Ey[w] + Ey[w_tp] );                              \
  dcomplex hx = -0.25*( Hx[w] + Hx[w_ft] + Hx[w_tp] + Hx[field_front(w_tp)]); \
  dcomplex hy =  0.25*( Hy[w] + Hy[w_rt] + Hy[w_ft] + Hy[field_front(w_rt)]);

//上下の面に必要なE,Hを求める
#define EH_IN_XZ(w, ex, ez, hx, hz)                                      \
  int w_rt = field_right(w);                                            \
  int w_tp = field_top(w);                                              \
  int w_ft = field_front(w);                                            \
  dcomplex ex =  0.5*( Ex[w] + Ex[w_rt] );                              \
  dcomplex ez = -0.5*( Ez[w] + Ez[w_ft] );                              \
  dcomplex hx = -0.25*( Hx[w] + Hx[w_ft] + Hx[w_tp] + Hx[field_front(w_tp)]); \
  dcomplex hz =  0.25*( Hz[w] + Hz[w_rt] + Hz[w_tp] + Hz[field_top(w_rt)]);

//左右の面に必要なE,Hを求める
#define EH_IN_YZ(w, ey,ez,hy,hz)            \
  int w_rt = field_right(w);               \
  int w_tp = field_top(w);                 \
  int w_ft = field_front(w);               \
  dcomplex ey =  0.5*( Ey[w] + Ey[w_tp] ); \
  dcomplex ez = -0.5*( Ez[w] + Ez[w_ft] ); \
  dcomplex hy =  0.25*( Hy[w] + Hy[w_rt] + Hy[w_ft] + Hy[field_front(w_rt)]); \
  dcomplex hz = -0.25*( Hz[w] + Hz[w_rt] + Hz[w_tp] + Hz[field_top(w_rt)]);


void ntff3D_Frequency( dcomplex *Ex, dcomplex *Ey,dcomplex *Ez,
                       dcomplex *Hx, dcomplex *Hy, dcomplex *Hz)
{
  NTFFInfo nInfo = field_getNTFFInfo();
  double cx = nInfo.cx;
  double cy = nInfo.cy;
  double cz = nInfo.cz;

  double R = 1.0e6;
  double k_s = field_getK();
  double w_s = field_getOmega();
  dcomplex Coeffician = I * w_s / (4*M_PI*R+C_0_S) * cexp(-I*w_s*R/C_0_S);
  int tp = nInfo.top;     //上面
  int bm = nInfo.bottom;  //下面
  int rt = nInfo.right;   //右
  int lt = nInfo.left;	  //左
  int ft = nInfo.front;   //全面
  int bk = nInfo.back;    //背面

  double Z0 = field_getZ_0_S();

  //phi   : z軸回転
  //theta : y(or x)軸回転 TODO
  //theta=0 => +x(右), theta=90=>+z(前)
  //phi = 90 => +y(上)
  double theta_rad = 0;
  for(int phi=0; phi<360; phi++) {
    double phi_rad = phi*M_PI/180.0;
    double r1x = cos(theta_rad)*cos(phi_rad), r1y = sin(phi_rad), r1z = cos(phi_rad)*sin(theta_rad);
    dcomplex Nx = 0, Ny = 0, Nz = 0;
    dcomplex Lx = 0, Ly = 0, Lz = 0;
    
    //前の面 n=(0, 0, 1)
    // (N =) J = n × H = ( -hy, hx, 0)
    // (L =) M = E × n = ( ey ,-ex,  0)
    for ( int i=lt; i<rt; i++ )
      for( int j=bm; j<tp; j++) {
        COEFF_AND_INDICES(i, j, ft);// coefやw を求める
        EH_IN_XY(w, ex, ey, hx, hy);   //ex, ey, hx, hyを補完して求める
        Lx += ey*coef;
        Ly -= ex*coef;
        
        Nx -= hy*coef;
        Ny += hx*coef;
      }

    //後の面 n=(0, 0, -1)
    // (N =) J = n × H = ( hy,-hx, 0)
    // (L =) M = E × n = (-ey, ex,  0)
    for ( int i=lt; i<rt; i++ )
      for( int j=bm; j<tp; j++) {
        COEFF_AND_INDICES(i, j, bk); //ファイルのトップ参照
        EH_IN_XY(w, ex, ey, hx, hy);   //ex, ey, hx, hyを補完して求める
        
        Lx -= ey*coef;
        Ly += ex*coef;
        Nx += hy*coef;
        Ny -= hx*coef;
      }

    //右の面 n=(1,0,0)
    // (N =) J = n × H = (  0, -hz, hy)
    // (L =) M = E × n = (  0,  ez,-ey)
    for(int j=bm; j<tp; j++)
      for(int k=bk; k<ft; k++){
        COEFF_AND_INDICES(rt, j, k);
        EH_IN_YZ(w, ey,ez, hy,hz);
        Ny -= hz*coef;
        Nz += hy*coef;
        Ly += ez*coef;
        Lz -= ey*coef;
      }

    //左の面 n=(-1, 0, 0)
    // (N =) J = n × H = (  0, hz,-hy)
    // (L =) M = E × n = (  0,-ez, ey)
    for(int j=bm; j<tp; j++)
      for(int k=bk; k<ft; k++){
        COEFF_AND_INDICES(lt, j, k);
        EH_IN_YZ(w, ey,ez, hy,hz);
        
        Ny += hz*coef;
        Nz -= hy*coef;
        Ly -= ez*coef;
        Lz += ey*coef;
      }

    //上の面 n=(0,1,0)
    // (N =) J = n × H = ( hz, 0, -hx)
    // (L =) M = E × n = (-ez, 0,  ex)
    for(int i=lt; i<rt; i++)
      for (int k=bk; k<ft; k++ )
      {
        COEFF_AND_INDICES(i, tp, k);
        EH_IN_XZ(w, ex, ez, hx, hz);

        Nx += hz*coef;
        Nz -= hx*coef;
        Lx -= ez*coef;
        Lz += ex*coef;
    }

    //下の面 n=(0,-1,0)
    // (N =) J = n × H = (-hz, 0, hx)
    // (L =) M = E × n = ( ez, 0,-ex)
    for(int i=lt; i<rt; i++)
      for (int k=bk; k<ft; k++ )
    {
      COEFF_AND_INDICES(i, bm, k);
      EH_IN_XZ(w, ex, ez, hx, hz);

      Nx -= hz*coef;
      Nz += hx*coef;
      Lx += ez*coef;
      Lz -= ex*coef;
    }

    double sx = cos(theta_rad)*cos(phi_rad);
    double sy = cos(theta_rad)*sin(phi_rad);
    double sz = -cos(theta_rad); //宇野先生の本では -sin(theta)になってる
    double px = -sin(phi_rad);
    double py = cos(phi_rad);

    dcomplex Nth = sx*Nx + sy*Ny + sz*Nz;
    dcomplex Nph = px*Nx + py*Ny;//-Nx*sin(phi_rad) + Ny*cos(phi_rad);
    dcomplex Lth = sx*Lx + sy*Ly + sz*Lz;
    dcomplex Lph = px*Lx + py*Ly;

    dcomplex Eth = Coeffician*(-Z0*Nth - Lph);
    dcomplex Eph = Coeffician*(-Z0*Nph + Lth);
  }
}


// 係数と最低限必要なインデックスをまとめて求める為のマクロ
#define SUB_TIME_SHIFT_AND_INDICES(i, j, k)                             \
  double r2x = i+0.5-cx, r2y = j+0.5-cy, r2z = k+0.5-cz;                \
  double dot_per_c = r1x_per_c*r2x + r1y_per_c*r2y + r1z_per_c*r2z;     \
  double timeShift = -dot_per_c + nInfo.RFperC;                         \

// eとU[stp] もしくは hとW[stp]を渡す
//UW_ang = Ux[stp],Uy[stp],Wz[stp] の事. array[360][num]を一次元配列で表しており, その角度における配列を引数にとる
static inline void calc(double time_plus_timeShift, dcomplex eh,  dcomplex *UW_ang){  
  int m = floor(time_plus_timeShift+0.5);
  double a = (0.5 + time_plus_timeShift - m);
  double b = 1.0-a;
  double ab = a-b;
  UW_ang[m-1] += eh*b;
  UW_ang[m]   += eh*ab;
  UW_ang[m+1] -= eh*a;
}

void ntff3D_SubTimeCalc(dcomplex *Ex,dcomplex *Ey,dcomplex *Ez,
                     dcomplex *Hx,dcomplex *Hy,dcomplex *Hz,
                     dcomplex *Ux,dcomplex *Uy,dcomplex *Uz,
                     dcomplex *Wx,dcomplex *Wy,dcomplex *Wz)
{
  double timeE = field_getTime() - 1;   //t - Δt
  double timeH = field_getTime() - 0.5; //t - Δt/2  
  NTFFInfo nInfo = field_getNTFFInfo();
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  double cx = nInfo.cx;
  double cy = nInfo.cy;
  double cz = nInfo.cz;

//  double R = 1.0e6;
//  double w_s = field_getOmega();
//  dcomplex Coeffician = I * w_s / (4*M_PI*R+C_0_S) * cexp(-I*w_s*R/C_0_S);
  
  int tp = nInfo.top;     //上面
  int bm = nInfo.bottom;  //下面
  int rt = nInfo.right;   //右
  int lt = nInfo.left;	  //左
  int ft = nInfo.front;   //全面
  int bk = nInfo.back;    //背面

  int subIndex_tp = tp - subInfo_s.OFFSET_Y;
  int subIndex_bm = bm - subInfo_s.OFFSET_Y;
  int subIndex_rt = rt - subInfo_s.OFFSET_X;
  int subIndex_lt = lt - subInfo_s.OFFSET_X;
  int subIndex_ft = ft - subInfo_s.OFFSET_Z;
  int subIndex_bk = bk - subInfo_s.OFFSET_Z;
  /*
  int sub_tp = min(subInfo_s.SUB_N_PY, max( 0, tp - subInfo_s.OFFSET_Y) );
  int sub_bm = min(subInfo_s.SUB_N_PY, max( 0, bm - subInfo_s.OFFSET_Y) );
  int sub_rt = min(subInfo_s.SUB_N_PX, max( 0, rt - subInfo_s.OFFSET_X) );
  int sub_lt = min(subInfo_s.SUB_N_PX, max( 0, lt - subInfo_s.OFFSET_X) );
  int sub_ft = min(subInfo_s.SUB_N_PZ, max( 0, ft - subInfo_s.OFFSET_Z) );
  int sub_bk = min(subInfo_s.SUB_N_PZ, max( 0, bk - subInfo_s.OFFSET_Z) );
  */
  // 小領域は上の積分面上にある
//  bool IN_TP =  (1 <= sub_tp != && sub_tp <= subInfo_s.N_PX-2) && sub_rt != sub_lt && sub_ft != sub_bk; // 
  
  int index_ang = 0;  //角度angの0番目のインデックス
  const double theta_rad = 0;
  const double ToRad = M_PI/180.0;
  for(int phi=0; phi < 360; phi++, index_ang+=nInfo.arraySize)
  {
    double phi_rad = phi*ToRad;
    double r1x_per_c = cos(theta_rad)*cos(phi_rad)/C_0_S;
    double r1y_per_c = sin(phi_rad)/C_0_S;
    double r1z_per_c = cos(phi_rad)*sin(theta_rad)/C_0_S;
    
    //ang°の位置にシフトしたポジション, こうすれば Ux_ang[i]でその角度のi番目にアクセスできる.
    dcomplex *Ux_ang = &Ux[index_ang];
    dcomplex *Uy_ang = &Uy[index_ang];
    dcomplex *Uz_ang = &Uz[index_ang];
    dcomplex *Wx_ang = &Wx[index_ang];
    dcomplex *Wy_ang = &Wy[index_ang];
    dcomplex *Wz_ang = &Wz[index_ang];

    //前の面 n=(0, 0, 1)
    // (N =) J = n × H = (-hy, hx, 0)
    // (L =) M = E × n = ( ey,-ex,  0)
    for ( int i=lt; i<rt; i++ )
      for( int j=bm; j<tp; j++) {
        SUB_TIME_SHIFT_AND_INDICES(i , j, ft);
        int w = field_subIndex(i-subInfo_s.OFFSET_X, j-subInfo_s.OFFSET_Y, ft-subInfo_s.OFFSET_Z);
        EH_IN_XY(w, ex, ey, hx, hy);
      
        calc(timeE+timeShift,-ex, Uy_ang);
        calc(timeE+timeShift, ey, Ux_ang);
        calc(timeH+timeShift, hx, Wy_ang);
        calc(timeH+timeShift,-hy, Wx_ang);
      }

    //前の面 n=(0, 0, 1)
    // (W =) J = n × H = (-hx, hy, 0)
    // (U =) M = E × n = ( ey,-ex,  0)
    for ( int i=lt; i<rt; i++ )
      for( int j=bm; j<tp; j++) {        
        SUB_TIME_SHIFT_AND_INDICES(i, j, bk);
        int w = field_subIndex(i-subInfo_s.OFFSET_X, j-subInfo_s.OFFSET_Y, bk-subInfo_s.OFFSET_Z);
        
        EH_IN_XY(w, ex, ey, hx, hy);
        calc(timeE+timeShift, ex, Uy_ang);
        calc(timeE+timeShift,-ey, Ux_ang);
        calc(timeH+timeShift,-hx, Wy_ang);
        calc(timeH+timeShift, hy, Wx_ang);
      }

    
    //右の面 n=(1,0,0)
    // (W =) J = n × H = (  0, -hz, hy)
    // (U =) M = E × n = (  0,  ez,-ey)
    for(int j=bm; j<tp; j++)
      for(int k=bk; k<ft; k++){
        SUB_TIME_SHIFT_AND_INDICES(rt, j, k);
        int w = field_subIndex(rt-subInfo_s.OFFSET_X, j-subInfo_s.OFFSET_Y, k-subInfo_s.OFFSET_Z);
        
        EH_IN_YZ(w, ey,ez,hy,hz);        
        calc(timeE+timeShift, ez, Uy_ang);
        calc(timeE+timeShift,-ey, Uz_ang);
        calc(timeH+timeShift,-hz, Wy_ang);
        calc(timeH+timeShift, hy, Wz_ang);
      }

    //左の面 n=(-1, 0, 0)
    // (W =) J = n × H = (  0, hz,-hy)
    // (U =) M = E × n = (  0,-ey, ez)
    for(int j=bm; j<tp; j++)
      for(int k=bk; k<ft; k++){
        SUB_TIME_SHIFT_AND_INDICES(lt, j, k);
        int w = field_subIndex(lt-subInfo_s.OFFSET_X, j-subInfo_s.OFFSET_Y, k-subInfo_s.OFFSET_Z);
        
        EH_IN_YZ(w, ey,ez, hy,hz);
        calc(timeE+timeShift,-ez, Uy_ang);
        calc(timeE+timeShift, ey, Uz_ang);
        calc(timeH+timeShift,-hz, Wy_ang);
        calc(timeH+timeShift, hy, Wz_ang);
      }

    //上の面 n=(0,1,0)
    // (W =) J = n × H = ( hz, 0, -hx)
    // (U =) M = E × n = (-ez, 0,  ex)
    for(int i=lt; i<rt; i++)
      for (int k=bk; k<ft; k++ )
      {
        SUB_TIME_SHIFT_AND_INDICES(i, tp, k);
        int w = field_subIndex(i-subInfo_s.OFFSET_X, tp-subInfo_s.OFFSET_Y, k-subInfo_s.OFFSET_Z);
        
        EH_IN_XZ(w, ex, ez, hx, hz);
        calc(timeE+timeShift,-ez, Ux_ang);
        calc(timeE+timeShift, ex, Uz_ang);
        calc(timeH+timeShift, hz, Wx_ang);
        calc(timeH+timeShift,-hx, Wz_ang);
      }

    //下の面 n=(0,-1,0)
    // (W =) J = n × H = (-hz, 0, hx)
    // (U =) M = E × n = ( ez, 0,-ex)
    for(int i=lt; i<rt; i++)
      for (int k=bk; k<ft; k++ )
      {
        SUB_TIME_SHIFT_AND_INDICES(i, bk, k);
        int w = field_subIndex(i-subInfo_s.OFFSET_X, bm-subInfo_s.OFFSET_Y, k-subInfo_s.OFFSET_Z);
        
        EH_IN_XZ(w, ex, ez, hx, hz);
        calc(timeE+timeShift, ez, Ux_ang);
        calc(timeE+timeShift,-ex, Uz_ang);
        calc(timeH+timeShift,-hz, Wx_ang);
        calc(timeH+timeShift, hx, Wz_ang);
      }
  }
}


// 以下 バックアップ用

// 周波数NTFF用
//前の面
/*
  double r2x = i+0.5-cx, r2y = j+0.5-cy, r2z = ft+0.5-cz;
  double dot = r1x*r2x + r1y*r2y + r1z*r2z;  //内積
  dcomplex coef = cexp( I*k_s*dot);
  int w = field_index(i, j, ft);
  int w_rt = field_right(w);
  int w_tp = field_top(w);
  int w_ft = field_front(w);
  int w_fttp = field_front(w_tp);
  int w_ftrt = field_front(w_rt);

  dcomplex ex = -0.5*( Ex[w] + Ex[w_rt] );
  dcomplex ey =  0.5*( Ey[w] + Ey[w_tp] );
  dcomplex hx = -0.25*( Hx[w] + Hx[w_ft] + Hx[w_tp] + Hx[w_fttp]);
  dcomplex hy =  0.25*( Hy[w] + Hy[w_rt] + Hy[w_ft] + Hy[w_ftrt]);
*/

//後ろの面
/*
  double r2x = i+0.5-cx;  //セルの中心を積分路にする
  double r2y = j+0.5-cy;
  double r2z = bk+0.5-cz;
  double dot = r1x*r2x + r1y*r2y + r1z*r2z;  //内積
  dcomplex coef = cexp( I*k_s*dot);
  int w = field_index(i, j, bk);
  int w_rt = field_right(w);
  int w_tp = field_top(w);
  int w_ft = field_front(w);  
  int w_fttp = field_front(w_tp);
  int w_ftrt = field_front(w_rt);

  dcomplex ex =  0.5*( Ex[w] + Ex[w_rt] );
  dcomplex ey = -0.5*( Ey[w] + Ey[w_tp] );
  dcomplex hx =  0.25*( Hx[w] + Hx[w_ft] + Hx[w_tp] + Hx[w_fttp]);
  dcomplex hy = -0.25*( Hy[w] + Hy[w_rt] + Hy[w_ft] + Hy[w_ftrt]);
*/

//右の面
/*
  double r2x = rt+0.5-cx;
  double r2y = j+0.5-cy;
  double r2z = k+0.5-cz;
  double dot = r1x*r2x + r1y*r2y + r1z*r2z;  //内積
  dcomplex coef = cexp( I*k_s*dot);        
  int w = field_index(rt, j, k);
  int w_rt = field_right(w);
  int w_tp = field_top(w);
  int w_ft = field_front(w);          

  int w_tprt = field_top(w_rt);
  int w_ftrt = field_front(w_rt);
  dcomplex ey =  0.5*( Ey[w] + Ey[w_tp] );
  dcomplex ez = -0.5*( Ez[w] + Ez[w_ft] );
  dcomplex hy =  0.25*( Hy[w] + Hy[w_rt] + Hy[w_ft] + Hy[w_ftrt]);        
  dcomplex hz = -0.25*( Hz[w] + Hz[w_rt] + Hz[w_tp] + Hz[w_tprt]);
*/

//左の面
/*
  double r2x = lt+0.5-cx;
  double r2y = j+0.5-cy;
  double r2z = k+0.5-cz;
  double dot = r1x*r2x + r1y*r2y + r1z*r2z;  //内積
  dcomplex coef = cexp( I*k_s*dot);        
  int w = field_index(lt, j, k);
  int w_rt = field_right(w);
  int w_tp = field_top(w);
  int w_ft = field_front(w);
  int w_tprt = field_top(w_rt);
  int w_ftrt = field_front(w_rt);
  dcomplex ey = -0.5*( Ey[w] + Ey[w_tp] );
  dcomplex ez =  0.5*( Ez[w] + Ez[w_ft] );
  dcomplex hy = -0.25*( Hy[w] + Hy[w_rt] + Hy[w_ft] + Hy[w_ftrt]);        
  dcomplex hz =  0.25*( Hz[w] + Hz[w_rt] + Hz[w_tp] + Hz[w_tprt]);
*/

//上の面
/*
  double r2x = i+0.5-cx;
  double r2y = tp+0.5-cy;
  double r2z = k+0.5-cz;
  double dot = r1x*r2x + r1y*r2y + r1z*r2z;  //内積
  dcomplex coef = cexp( I*k_s*dot);
  int w = field_index(i, tp, k);
  int w_rt = field_right(w);
  int w_tp = field_top(w);
  int w_ft = field_front(w);
  int w_tprt = field_top(w_rt);
  int w_fttp = field_front(w_tp);

  dcomplex ex =  0.5*( Ex[w] + Ex[w_rt] );
  dcomplex ez = -0.5*( Ez[w] + Ez[w_ft] );
  dcomplex hx = -0.25*( Hx[w] + Hx[w_ft] + Hx[w_tp] + Hx[w_fttp]);
  dcomplex hz =  0.25*( Hz[w] + Hz[w_rt] + Hz[w_tp] + Hz[w_tprt]);
*/

//下の面
/*
  double r2x =  i+0.5-cx;
  double r2y = bm+0.5-cy;
  double r2z =  k+0.5-cz;
  double dot = r1x*r2x + r1y*r2y + r1z*r2z;  //内積
  dcomplex coef = cexp( I*k_s*dot);
  int w = field_index(i, bm, k);
  int w_rt = field_right(w);
  int w_tp = field_top(w);
  int w_ft = field_front(w);
  int w_tprt = field_top(w_rt);
  int w_fttp = field_front(w_tp);

  dcomplex ex =  0.5*( Ex[w] + Ex[w_rt] );
  dcomplex ez = -0.5*( Ez[w] + Ez[w_ft] );
  dcomplex hx = -0.25*( Hx[w] + Hx[w_ft] + Hx[w_tp] + Hx[w_fttp]);
  dcomplex hz =  0.25*( Hz[w] + Hz[w_rt] + Hz[w_tp] + Hz[w_tprt]);

*/
