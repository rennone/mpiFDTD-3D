#include <math.h>
#include <stdlib.h>
#include "ntff3D.h"
#include "field.h"
#include "function.h"
#include "myComplex.h"

static dcomplex *Ux,*Uy,*Uz,*Wx,*Wy,*Wz;
static double R; //観測点までの距離

//sub領域に置ける, 積分路
static int sub_ylt, sub_yrt, sub_yft, sub_ybk; //tp,bm面に置ける左右前後の範囲(直方体なので, tp面とbm面で違いは無い)
static int sub_xtp, sub_xbm, sub_xft, sub_xbk;
static int sub_zlt, sub_zrt, sub_ztp, sub_zbm;

static int sub_tp, sub_bm, sub_rt, sub_lt, sub_ft, sub_bk;
static bool IN_TP, IN_BM, IN_LT, IN_RT, IN_FT, IN_BK;

void ntff3D_Init()
{
  NTFFInfo nInfo = field_getNTFFInfo();

  
  Ux = newDComplex(nInfo.arraySize * 360);
  Uy = newDComplex(nInfo.arraySize * 360);
  Uz = newDComplex(nInfo.arraySize * 360);
  Wx = newDComplex(nInfo.arraySize * 360);
  Wy = newDComplex(nInfo.arraySize * 360);
  Wz = newDComplex(nInfo.arraySize * 360);

  R = 1.0e6 * field_toCellUnit(500);//* field_getLambda_S();

  double cx = nInfo.cx;
  double cy = nInfo.cy;
  double cz = nInfo.cz;
  
  int tp = nInfo.top;    int bm = nInfo.bottom;  //上下
  int rt = nInfo.right;  int lt = nInfo.left;	 //左右
  int ft = nInfo.front;  int bk = nInfo.back;    //前後

  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  
  sub_tp = tp - subInfo_s.OFFSET_Y;  sub_bm = bm - subInfo_s.OFFSET_Y;
  sub_rt = rt - subInfo_s.OFFSET_X;  sub_lt = lt - subInfo_s.OFFSET_X;
  sub_ft = ft - subInfo_s.OFFSET_Z;  sub_bk = bk - subInfo_s.OFFSET_Z;

  //以下どれかでも満たせば積分路上に無い
  bool outX = sub_rt <= 0 || sub_lt >= subInfo_s.SUB_N_PX-1; //rtより右, もしくはltより左の小領域
  bool outY = sub_tp <= 0 || sub_bm >= subInfo_s.SUB_N_PY-1; //tpより上, もしくはbmより下の小領域
  bool outZ = sub_ft <= 0 || sub_bk >= subInfo_s.SUB_N_PZ-1; //ftより前, もしくはbkより後ろの小領域
  
  // 小領域内にどの積分面が存在するか
  IN_TP = (0 < sub_tp && sub_tp < subInfo_s.SUB_N_PY-1) && !outX && !outZ;
  IN_BM = (0 < sub_bm && sub_bm < subInfo_s.SUB_N_PY-1) && !outX && !outZ;
  IN_RT = (0 < sub_rt && sub_rt < subInfo_s.SUB_N_PX-1) && !outY && !outZ;
  IN_LT = (0 < sub_lt && sub_lt < subInfo_s.SUB_N_PX-1) && !outY && !outZ;
  IN_FT = (0 < sub_ft && sub_ft < subInfo_s.SUB_N_PZ-1) && !outX && !outY;
  IN_BK = (0 < sub_bk && sub_bk < subInfo_s.SUB_N_PZ-1) && !outX && !outY;

  sub_ylt=-1, sub_yrt=-2, sub_ybk=-1, sub_yft=-2 ;
  if(IN_TP || IN_BM)
  {
    sub_yrt = min(subInfo_s.SUB_N_PX-2, max( 1, sub_rt) );
    sub_ylt = min(subInfo_s.SUB_N_PX-2, max( 1, sub_lt) );
    sub_yft = min(subInfo_s.SUB_N_PZ-2, max( 1, sub_ft) );
    sub_ybk = min(subInfo_s.SUB_N_PZ-2, max( 1, sub_bk) );
  }

  sub_xtp=-2, sub_xbm=-1, sub_xbk=-1, sub_xft=-2;
  if(IN_RT || IN_LT)
  {
    sub_xbm = min(subInfo_s.SUB_N_PY-2, max( 1, sub_bm+1) );  //bm,tpですでに計算しているため, ひとつずれる
    sub_xtp = min(subInfo_s.SUB_N_PY-2, max( 1, sub_tp-1) );  //
    sub_xft = min(subInfo_s.SUB_N_PZ-2, max( 1, sub_ft) );
    sub_xbk = min(subInfo_s.SUB_N_PZ-2, max( 1, sub_bk) );
  }

  sub_ztp=-2, sub_zbm=-1, sub_zlt=-1, sub_zrt=-2;
  if(IN_FT || IN_BK)
  {
    sub_zbm = min(subInfo_s.SUB_N_PY-2, max( 1, sub_bm+1) );  //bm,tp面ですでに計算しているため,一つずれる
    sub_ztp = min(subInfo_s.SUB_N_PY-2, max( 1, sub_tp-1) );
    sub_zrt = min(subInfo_s.SUB_N_PX-2, max( 1, sub_rt-1) );  //lt,rtも同様
    sub_zlt = min(subInfo_s.SUB_N_PX-2, max( 1, sub_lt+1) );
  }

  printf("Rank=%d\ninTp=%d, inBm=%d, inLt=%d, inRt=%d, inFt=%d, inBk=%d\n offset(%d,%d,%d)\n",subInfo_s.Rank,IN_TP, IN_BM, IN_RT,IN_LT,IN_FT, IN_BK, subInfo_s.OFFSET_X, subInfo_s.OFFSET_Y, subInfo_s.OFFSET_Z);

}

//必ず COEFF_AND_INDICESの後に記述する
//前後の面に必要なE,Hを求める
#define EH_IN_XY(w, ex, ey, hx, hy)                                      \
  int w_rt = field_right(w);                                            \
  int w_tp = field_top(w);                                              \
  int w_ft = field_front(w);                                            \
  dcomplex ex =  0.5*( Ex[w] + Ex[w_rt] );                              \
  dcomplex ey =  0.5*( Ey[w] + Ey[w_tp] );                              \
  dcomplex hx =  0.25*( Hx[w] + Hx[w_ft] + Hx[w_tp] + Hx[field_front(w_tp)]); \
  dcomplex hy =  0.25*( Hy[w] + Hy[w_rt] + Hy[w_ft] + Hy[field_front(w_rt)]);

//上下の面に必要なE,Hを求める
#define EH_IN_XZ(w, ex, ez, hx, hz)                                      \
  int w_rt = field_right(w);                                            \
  int w_tp = field_top(w);                                              \
  int w_ft = field_front(w);                                            \
  dcomplex ex =  0.5*( Ex[w] + Ex[w_rt] );                              \
  dcomplex ez =  0.5*( Ez[w] + Ez[w_ft] );                              \
  dcomplex hx =  0.25*( Hx[w] + Hx[w_ft] + Hx[w_tp] + Hx[field_front(w_tp)]); \
  dcomplex hz =  0.25*( Hz[w] + Hz[w_rt] + Hz[w_tp] + Hz[field_top(w_rt)]);

//左右の面に必要なE,Hを求める
#define EH_IN_YZ(w, ey,ez,hy,hz)            \
  int w_rt = field_right(w);               \
  int w_tp = field_top(w);                 \
  int w_ft = field_front(w);               \
  dcomplex ey =  0.5*( Ey[w] + Ey[w_tp] ); \
  dcomplex ez =  0.5*( Ez[w] + Ez[w_ft] ); \
  dcomplex hy =  0.25*( Hy[w] + Hy[w_rt] + Hy[w_ft] + Hy[field_front(w_rt)]); \
  dcomplex hz =  0.25*( Hz[w] + Hz[w_rt] + Hz[w_tp] + Hz[field_top(w_rt)]);

#define CALC_COEFF(i, j, k)                      \
  double r2x = i+0.5-cx, r2y = j+0.5-cy, r2z = k+0.5-cz; \
  double dot = r1x*r2x + r1y*r2y + r1z*r2z;               \
  dcomplex coef = cexp( I*k_s*dot);                     \



static void frequencyNTFF(dcomplex *Ex, dcomplex *Ey,dcomplex *Ez,
                          dcomplex *Hx, dcomplex *Hy, dcomplex *Hz,
                          dcomplex *Eth, dcomplex *Eph, int theta, int phi)
{
  double k_s = field_getK_S();
  dcomplex Coeffician = I * k_s / (4*M_PI*R) * cexp(-I*k_s*R); //
  NTFFInfo nInfo = field_getNTFFInfo();
  double cx = nInfo.cx;  double cy = nInfo.cy;  double cz = nInfo.cz;
  int tp = nInfo.top;    int bm = nInfo.bottom; //上下
  int rt = nInfo.right;  int lt = nInfo.left;   //左右
  int ft = nInfo.front;  int bk = nInfo.back;   //前後
  double Z0 = field_getZ_0_S();
  double ToRad = M_PI/180.0;
  
  double theta_rad = theta * ToRad;
  double phi_rad   =   phi * ToRad;
  double r1x = sin(theta_rad)*cos(phi_rad);
  double r1y = sin(theta_rad)*sin(phi_rad);
  double r1z = cos(theta_rad);

  dcomplex Nx = 0, Ny = 0, Nz = 0;
  dcomplex Lx = 0, Ly = 0, Lz = 0;

  //上下は, 蓋をするように全範囲を網羅
  //上の面 n=(0,1,0)
  // (N =) J = n × H = ( hz, 0, -hx)
// (L =) M = E × n = (-ez, 0,  ex)
  for(int i=lt; i<=rt; i++)
    for (int k=bk; k<=ft; k++ )
    {      
      CALC_COEFF(i, tp, k);
      int w = field_index(i, tp, k);
      EH_IN_XZ(w, ex, ez, hx, hz);
      Nx += hz*coef;      Nz -= hx*coef;
      Lx -= ez*coef;      Lz += ex*coef;
    }
  
  //下の面 n=(0,-1,0)
  // (N =) J = n × H = (-hz, 0, hx)
  // (L =) M = E × n = ( ez, 0,-ex)
  for(int i=lt; i<=rt; i++)
    for (int k=bk; k<=ft; k++ )
    {
      CALC_COEFF(i, bm, k);
      int w = field_index(i, bm, k);
      EH_IN_XZ(w, ex, ez, hx,hz);
      Nx -= hz*coef;      Nz += hx*coef;
      Lx += ez*coef;      Lz -= ex*coef;
    }

  //前の面 n=(0, 0, 1)
  // (N =) J = n × H = ( -hy, hx, 0)
  // (L =) M = E × n = ( ey ,-ex,  0)
  for ( int i=lt; i<rt; i++ )
    for( int j=bm+1; j<tp; j++) {
      CALC_COEFF(i, j, ft);
      int w = field_index(i, j, ft);
      EH_IN_XY(w, ex, ey, hx, hy);//ex, ey, hx, hyを補完して求める
      Nx -= hy*coef;      Ny += hx*coef;
      Lx += ey*coef;      Ly -= ex*coef;
    }

  //後の面 n=(0, 0, -1)
  // (N =) J = n × H = ( hy,-hx, 0)
  // (L =) M = E × n = (-ey, ex,  0)
  for ( int i=lt; i<rt; i++ )
    for( int j=bm+1; j<tp; j++) {
      CALC_COEFF(i, j, bk); //このファイルの最初に書いてある
      int w = field_index(i, j, bk);
      EH_IN_XY(w, ex, ey, hx, hy);//ex, ey, hx, hyを補完して求める
      Nx += hy*coef;      Ny -= hx*coef;
      Lx -= ey*coef;      Ly += ex*coef;      
    }

  //右の面 n=(1,0,0)
  // (N =) J = n × H = (  0, -hz, hy)
  // (L =) M = E × n = (  0,  ez,-ey)
  for(int j=bm+1; j<tp; j++)
    for(int k=bk; k<ft; k++){
      CALC_COEFF(rt, j, k);
      int w = field_index(rt, j, k);
      EH_IN_YZ(w, ey,ez, hy,hz);
      Ny -= hz*coef;      Nz += hy*coef;
      Ly += ez*coef;      Lz -= ey*coef;
    }

  //左の面 n=(-1, 0, 0)
  // (N =) J = n × H = (  0, hz,-hy)
  // (L =) M = E × n = (  0,-ez, ey)
  for(int j=bm+1; j<tp; j++)
    for(int k=bk+1; k<ft; k++){
      CALC_COEFF(lt, j, k);
      int w = field_index(lt, j, k);
      EH_IN_YZ(w, ey,ez, hy,hz);
      Ny += hz*coef;      Nz -= hy*coef;
      Ly -= ez*coef;      Lz += ey*coef;
    }
      
  double sx = cos(theta_rad)*cos(phi_rad);
  double sy = cos(theta_rad)*sin(phi_rad);
  double sz = -sin(theta_rad);
  double px = -sin(phi_rad);
  double py = cos(phi_rad);

  dcomplex Nth = sx*Nx + sy*Ny + sz*Nz;
  dcomplex Nph = px*Nx + py*Ny;
  dcomplex Lth = sx*Lx + sy*Ly + sz*Lz;
  dcomplex Lph = px*Lx + py*Ly;

  *Eth = Coeffician*(Z0*Nth - Lph); //宇野先生の本では Nth, Nphに-がついていた.
  *Eph = Coeffician*(Z0*Nph + Lth);
}

void ntff3D_Frequency( dcomplex *Ex, dcomplex *Ey,dcomplex *Ez,
                       dcomplex *Hx, dcomplex *Hy, dcomplex *Hz)
{
  dcomplex Eth, Eph;
    
  //XY平面の遠方解を出力
  {
    FieldInfo fInfo = field_getFieldInfo();
    char buf[512];
    sprintf(buf, "%dnm_Eth_str.txt", fInfo.lambda_nm);
    FILE *fpEth = openFile(buf);
    sprintf(buf, "%dnm_Eph_str.txt", fInfo.lambda_nm);
    FILE *fpEph = openFile(buf);
    for(int theta=0; theta<360; theta++)
    {
      for(int phi=0; phi<360; phi++)
      {
        frequencyNTFF(Ex, Ey, Ez, Hx, Hy, Hz, &Eth, &Eph, theta, phi);
        fprintf(fpEth, "%.20lf ", cnorm(Eth));
        fprintf(fpEph, "%.20lf ", cnorm(Eph));
      }
      fprintf(fpEth, "\n");
      fprintf(fpEph, "\n");
    }
    fclose(fpEth);
    fclose(fpEph);
  }
}

// 係数と最低限必要なインデックスをまとめて求める為のマクロ
#define SUB_TIME_SHIFT(i, j, k)                             \
  double r2x = i+0.5-cx, r2y = j+0.5-cy, r2z = k+0.5-cz;                \
  double dot_per_c = r1x_per_c*r2x + r1y_per_c*r2y + r1z_per_c*r2z;     \
  double timeShift = -dot_per_c + nInfo.RFperC;                         \

//必ず COEFF_AND_INDICESの後に記述する
//前後の面に必要なE,Hを求める
#define SUB_EH_IN_XY(w, ex, ey, hx, hy)                                      \
  int w_rt = field_subRight(w);                                            \
  int w_tp = field_subTop(w);                                              \
  int w_ft = field_subFront(w);                                            \
  dcomplex ex = 0.5*( Ex[w] + Ex[w_rt] );                              \
  dcomplex ey = 0.5*( Ey[w] + Ey[w_tp] );                              \
  dcomplex hx = 0.25*( Hx[w] + Hx[w_ft] + Hx[w_tp] + Hx[field_front(w_tp)]); \
  dcomplex hy = 0.25*( Hy[w] + Hy[w_rt] + Hy[w_ft] + Hy[field_front(w_rt)]);

//上下の面に必要なE,Hを求める
#define SUB_EH_IN_XZ(w, ex, ez, hx, hz)                                      \
  int w_rt = field_subRight(w);                                            \
  int w_tp = field_subTop(w);                                              \
  int w_ft = field_subFront(w);                                            \
  dcomplex ex = 0.5*( Ex[w] + Ex[w_rt] );                              \
  dcomplex ez = 0.5*( Ez[w] + Ez[w_ft] );                              \
  dcomplex hx = 0.25*( Hx[w] + Hx[w_ft] + Hx[w_tp] + Hx[field_front(w_tp)]); \
  dcomplex hz = 0.25*( Hz[w] + Hz[w_rt] + Hz[w_tp] + Hz[field_top(w_rt)]);

//左右の面に必要なE,Hを求める
#define SUB_EH_IN_YZ(w, ey,ez,hy,hz)            \
  int w_rt = field_subRight(w);               \
  int w_tp = field_subTop(w);                 \
  int w_ft = field_subFront(w);               \
  dcomplex ey = 0.5*( Ey[w] + Ey[w_tp] ); \
  dcomplex ez = 0.5*( Ez[w] + Ez[w_ft] ); \
  dcomplex hy = 0.25*( Hy[w] + Hy[w_rt] + Hy[w_ft] + Hy[field_front(w_rt)]); \
  dcomplex hz = 0.25*( Hz[w] + Hz[w_rt] + Hz[w_tp] + Hz[field_top(w_rt)]);

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
                               dcomplex *Hx,dcomplex *Hy,dcomplex *Hz)
{
  double timeE = field_getTime() - 1;   //t - Δt
  double timeH = field_getTime() - 0.5; //t - Δt/2  
  NTFFInfo nInfo = field_getNTFFInfo();
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  double cx = nInfo.cx;
  double cy = nInfo.cy;
  double cz = nInfo.cz;
  
  int tp = nInfo.top;    int bm = nInfo.bottom;  //上下
  int rt = nInfo.right;  int lt = nInfo.left;	 //左右
  int ft = nInfo.front;  int bk = nInfo.back;    //前後
  
  int index_ang = 0;  //角度angの0番目のインデックス

  double phi_rad = M_PI/2.0;
  const double ToRad = M_PI/180.0;
  for(int theta=0; theta < 360; theta++, index_ang+=nInfo.arraySize)
  {
    double theta_rad = theta*ToRad;
    double r1x_per_c = sin(theta_rad)*cos(phi_rad)/C_0_S;
    double r1y_per_c = sin(theta_rad)*sin(phi_rad)/C_0_S;
    double r1z_per_c = cos(theta_rad)/C_0_S;
    
    //ang°の位置にシフトしたポジション, こうすれば Ux_ang[i]でその角度のi番目にアクセスできる.
    dcomplex *Ux_ang = &Ux[index_ang];    dcomplex *Uy_ang = &Uy[index_ang];
    dcomplex *Uz_ang = &Uz[index_ang];    dcomplex *Wx_ang = &Wx[index_ang];
    dcomplex *Wy_ang = &Wy[index_ang];    dcomplex *Wz_ang = &Wz[index_ang];

    //前の面 n=(0, 0, 1)
    // (W =) J = n × H = (-hy, hx, 0)
    // (U =) M = E × n = ( ey,-ex,  0)
    if( IN_FT )
    {
      for ( int i=sub_zlt; i<=sub_zrt; i++ )
        for( int j=sub_zbm; j<=sub_ztp; j++) {
          SUB_TIME_SHIFT(i+subInfo_s.OFFSET_X , j+subInfo_s.OFFSET_Y, ft);
          int w = field_subIndex(i, j, sub_ft);
          SUB_EH_IN_XY(w, ex, ey, hx, hy);
      
          calc(timeE+timeShift,-ex, Uy_ang);
          calc(timeE+timeShift, ey, Ux_ang);
          calc(timeH+timeShift, hx, Wy_ang);
          calc(timeH+timeShift,-hy, Wx_ang);
        }
    }
    
    //前の面 n=(0, 0, 1)
    // (W =) J = n × H = (-hx, hy, 0)
    // (U =) M = E × n = ( ey,-ex,  0)
    if( IN_BK )
    {
      for ( int i=sub_zlt; i<=sub_zrt; i++ )
        for( int j=sub_zbm; j<=sub_ztp; j++) {        
          SUB_TIME_SHIFT(i+subInfo_s.OFFSET_X, j+subInfo_s.OFFSET_Y, bk);
          int w = field_subIndex(i, j, sub_bk);
        
          SUB_EH_IN_XY(w, ex, ey, hx, hy);
          calc(timeE+timeShift, ex, Uy_ang);
          calc(timeE+timeShift,-ey, Ux_ang);
          calc(timeH+timeShift,-hx, Wy_ang);
          calc(timeH+timeShift, hy, Wx_ang);
        }
    }
    
    //右の面 n=(1,0,0)
    // (W =) J = n × H = (  0, -hz, hy)
    // (U =) M = E × n = (  0,  ez,-ey)
    if( IN_RT )
    {
      for(int j=sub_xbm; j<=sub_xtp; j++)
        for(int k=sub_xbk; k<=sub_xft; k++){
          SUB_TIME_SHIFT(sub_rt, j+subInfo_s.OFFSET_Y, k+subInfo_s.OFFSET_Z);
          int w = field_subIndex(sub_rt, j, k);
        
          SUB_EH_IN_YZ(w, ey,ez,hy,hz);        
          calc(timeE+timeShift, ez, Uy_ang);
          calc(timeE+timeShift,-ey, Uz_ang);
          calc(timeH+timeShift,-hz, Wy_ang);
          calc(timeH+timeShift, hy, Wz_ang);
        }
    }
    
    //左の面 n=(-1, 0, 0)
    // (W =) J = n × H = (  0, hz,-hy)
    // (U =) M = E × n = (  0,-ey, ez)
    if( IN_LT )
    {
      for(int j=sub_xbm; j<sub_xtp; j++)
        for(int k=sub_xbk; k<sub_xft; k++){
          SUB_TIME_SHIFT(sub_lt, j+subInfo_s.OFFSET_Y, k+subInfo_s.OFFSET_Z);
          int w = field_subIndex(sub_lt, j, k);
        
          SUB_EH_IN_YZ(w, ey,ez, hy,hz);
          calc(timeE+timeShift,-ez, Uy_ang);
          calc(timeE+timeShift, ey, Uz_ang);
          calc(timeH+timeShift,-hz, Wy_ang);
          calc(timeH+timeShift, hy, Wz_ang);
        }
    }
    
    //上の面 n=(0,1,0)
    // (W =) J = n × H = ( hz, 0, -hx)
    // (U =) M = E × n = (-ez, 0,  ex)
    if( IN_TP )
    {
      for(int i=sub_ylt; i<=sub_yrt; i++)
        for (int k=sub_ybk; k<=sub_yft; k++ )
        {
          SUB_TIME_SHIFT(i+subInfo_s.OFFSET_X, tp, k+subInfo_s.OFFSET_Z);
          int w = field_subIndex(i, sub_tp, k);
        
          SUB_EH_IN_XZ(w, ex, ez, hx, hz);
          calc(timeE+timeShift,-ez, Ux_ang);
          calc(timeE+timeShift, ex, Uz_ang);
          calc(timeH+timeShift, hz, Wx_ang);
          calc(timeH+timeShift,-hx, Wz_ang);
        }
    }
    
    //下の面 n=(0,-1,0)
    // (W =) J = n × H = (-hz, 0, hx)
    // (U =) M = E × n = ( ez, 0,-ex)
    if( IN_BM )
    {
      for(int i=sub_ylt; i<=sub_yrt; i++)
        for (int k=sub_ybk; k<=sub_yft; k++ )
        {
          SUB_TIME_SHIFT(i+subInfo_s.OFFSET_X, bk, k+subInfo_s.OFFSET_Z);
          int w = field_subIndex(i, sub_bm, k);
        
          EH_IN_XZ(w, ex, ez, hx, hz);
          calc(timeE+timeShift, ez, Ux_ang);
          calc(timeE+timeShift,-ex, Uz_ang);
          calc(timeH+timeShift,-hz, Wx_ang);
          calc(timeH+timeShift, hx, Wz_ang);
        }
    }
  }
}

static void ntff3D_TimeTranslate(dcomplex *Ux_ang, dcomplex *Uy_ang, dcomplex *Uz_ang,
                                 dcomplex *Wx_ang, dcomplex *Wy_ang, dcomplex *Wz_ang,
                                 dcomplex *Eth, dcomplex *Eph,
                                 int theta, int phi)
{
  const double complex coef = 1.0/(4*M_PI*C_0_S*R);
  const int maxTime = field_getMaxTime();

  double ToRad = M_PI/180.0;

  double theta_rad = theta * ToRad;
  double phi_rad   = phi*ToRad;
  double sx = cos(theta_rad)*cos(phi_rad);
  double sy = cos(theta_rad)*sin(phi_rad);
  double sz = -sin(theta_rad); //宇野先生の本では -sin(theta)になってる
  double px = -sin(phi_rad);
  double py = cos(phi_rad);
  double Z0 = field_getZ_0_S();
  for(int i=0; i < maxTime; i++)
  {
    double complex WTH = Wx_ang[i]*sx + Wy_ang[i]*sy + Wz_ang[i]*sz;
    double complex WPH = Wx_ang[i]*px + Wy_ang[i]*sy;
    double complex UTH = Ux_ang[i]*sx + Uy_ang[i]*sy + Uz_ang[i]*sz;
    double complex UPH = Ux_ang[i]*px + Uy_ang[i]*py;
    double complex ETH = coef*(-Z0*WTH-UPH);
    double complex EPH = coef*(-Z0*WPH+UTH);
      
    Eth[i] = ETH; //TODO : 物理単位に変換
    Eph[i] = EPH;
  }

}

static void unifyToRank0(dcomplex *p)
{
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  NTFFInfo nInfo = field_getNTFFInfo();
  int size = 360*nInfo.arraySize;

  if( subInfo_s.Rank == 0 )
  {
    MPI_Status status;
    dcomplex *tmp = newDComplex(size);
    for(int i=1; i<subInfo_s.Nproc; i++)
    {
      MPI_Recv(tmp, size, MPI_C_DOUBLE_COMPLEX, i, 0, MPI_COMM_WORLD, &status);
      for(int j=0; j<size; j++)
        p[j] += tmp[j];
    }
    free(tmp);
  }
  else {
    MPI_Send(p, size, MPI_C_DOUBLE_COMPLEX, 0, 0, MPI_COMM_WORLD);
  }
}

// データの書き出し.
void outputTimeDomainData(const char *fileName, dcomplex *data)
{
  char name[512];
  sprintf(name, "re_%s",fileName);  FILE *reFp = openFile(name);
  sprintf(name, "im_%s",fileName);  FILE *imFp = openFile(name);

  NTFFInfo nInfo = field_getNTFFInfo();
  const int maxTime = field_getMaxTime();
  for(int ang=0; ang<360; ang++)
  {
    int k= ang*nInfo.arraySize;
    for(int i=0; i < maxTime; i++)
    {
      fprintf(reFp,"%.20lf " , creal(data[k+i]));
      fprintf(imFp,"%.20lf " , cimag(data[k+i]));
    }
    fprintf(reFp,"\n");
    fprintf(imFp,"\n");
  }
  fclose(reFp);
  fclose(imFp);
}

//時間領域のEthの書き出し.
void ntff3D_TimeOutput()
{
  NTFFInfo nInfo = field_getNTFFInfo();
  dcomplex *Eth, *Eph;
  int size = 360*nInfo.arraySize;
  Eth = newDComplex(size);
  Eph = newDComplex(size);

  for(int theta=0; theta<360; theta++)
  {
    int ind = theta * nInfo.arraySize;
    ntff3D_TimeTranslate(&Ux[ind], &Uy[ind], &Uz[ind],
                         &Wx[ind], &Wy[ind], &Wz[ind],
                         &Eth[ind], &Eph[ind],
                         theta, 90);
  }

  //ランク0に集める
  unifyToRank0(Ux);  unifyToRank0(Uy);    unifyToRank0(Uz);
  unifyToRank0(Wx);  unifyToRank0(Wy);    unifyToRank0(Wz);
  unifyToRank0(Eth); unifyToRank0(Eph);
  SubFieldInfo_S sInfo_s = field_getSubFieldInfo_S();
  if(sInfo_s.Rank == 0)
  {

    outputTimeDomainData("Ux_zy.txt", Ux);
    outputTimeDomainData("Uy_zy.txt", Uy);
    outputTimeDomainData("Uz_zy.txt", Uz);
    outputTimeDomainData("Wx_zy.txt", Wx);
    outputTimeDomainData("Wy_zy.txt", Wy);
    outputTimeDomainData("Wz_zy.txt", Wz);
    outputTimeDomainData("Eth_time_zy.txt", Eth);
    outputTimeDomainData("Eph_time_zy.txt", Eph);
  }
  
  free(Eth);
  free(Eph);
}

//小領域に置ける,周波数NTFF
static inline void subFrequencyNTFF(dcomplex *Ex, dcomplex *Ey,dcomplex *Ez,
                             dcomplex *Hx, dcomplex *Hy, dcomplex *Hz,
                             dcomplex *Eth, dcomplex *Eph,
                                    double theta_rad, double phi_rad, dcomplex Coeffician)
{
  double k_s = field_getK_S();
//  dcomplex Coeffician =  I * k_s / (4*M_PI*R) * cexp(-I*k_s*R);
  NTFFInfo nInfo = field_getNTFFInfo();
  double cx = nInfo.cx;  double cy = nInfo.cy;  double cz = nInfo.cz;
  int tp = nInfo.top;    int bm = nInfo.bottom; //上下
  int rt = nInfo.right;  int lt = nInfo.left;   //左右
  int ft = nInfo.front;  int bk = nInfo.back;   //前後
  double Z0 = field_getZ_0_S();

  double r1x = sin(theta_rad)*cos(phi_rad);
  double r1y = sin(theta_rad)*sin(phi_rad);
  double r1z = cos(theta_rad);

  dcomplex Nx = 0, Ny = 0, Nz = 0;
  dcomplex Lx = 0, Ly = 0, Lz = 0;

  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  //上下は, 蓋をするように全範囲を網羅
  //上の面 n=(0,1,0)
  // (N =) J = n × H = ( hz, 0, -hx)
// (L =) M = E × n = (-ez, 0,  ex)
  if(IN_TP)
  {
    for(int i=sub_ylt; i<=sub_yrt; i++)
      for (int k=sub_ybk; k<=sub_yft; k++ )
      {      
        CALC_COEFF(i+subInfo_s.OFFSET_X, tp, k+subInfo_s.OFFSET_Z);
        int w = field_subIndex(i, sub_tp, k);
        SUB_EH_IN_XZ(w, ex, ez, hx, hz);
        Nx += hz*coef;      Nz -= hx*coef;
        Lx -= ez*coef;      Lz += ex*coef;
      }
  }
  
  //下の面 n=(0,-1,0)
  // (N =) J = n × H = (-hz, 0, hx)
  // (L =) M = E × n = ( ez, 0,-ex)
  if(IN_BM)
  {
    for(int i=sub_ylt; i<=sub_yrt; i++)
      for (int k=sub_ybk; k<=sub_yft; k++ )
      {
        CALC_COEFF(i+subInfo_s.OFFSET_X, bm, k+subInfo_s.OFFSET_Z);
        int w = field_subIndex(i, sub_bm, k);
        SUB_EH_IN_XZ(w, ex, ez, hx,hz);
        Nx -= hz*coef;      Nz += hx*coef;
        Lx += ez*coef;      Lz -= ex*coef;
      }
  }
  
  //前の面 n=(0, 0, 1)
  // (N =) J = n × H = ( -hy, hx, 0)
  // (L =) M = E × n = ( ey ,-ex,  0)
  if(IN_FT)
  {
    for ( int i=sub_zlt; i<=sub_zrt; i++ )
      for( int j=sub_zbm; j<=sub_ztp; j++) {
        CALC_COEFF(i+subInfo_s.OFFSET_X, j+subInfo_s.OFFSET_Y, ft);
        int w = field_subIndex(i, j, sub_ft);
        SUB_EH_IN_XY(w, ex, ey, hx, hy);//ex, ey, hx, hyを補完して求める
        Nx -= hy*coef;      Ny += hx*coef;
        Lx += ey*coef;      Ly -= ex*coef;
      }
  }
  
  //後の面 n=(0, 0, -1)
  // (N =) J = n × H = ( hy,-hx, 0)
  // (L =) M = E × n = (-ey, ex,  0)
  if(IN_BK)
  {
    for ( int i=sub_zlt; i<=sub_zrt; i++ )
      for( int j=sub_zbm; j<=sub_ztp; j++) {
        CALC_COEFF(i+subInfo_s.OFFSET_X, j+subInfo_s.OFFSET_Y, bk); //このファイルの最初に書いてある
        int w = field_subIndex(i, j, sub_bk);
        SUB_EH_IN_XY(w, ex, ey, hx, hy);//ex, ey, hx, hyを補完して求める
        Nx += hy*coef;      Ny -= hx*coef;
        Lx -= ey*coef;      Ly += ex*coef;      
      }
  }
  
  //右の面 n=(1,0,0)
  // (N =) J = n × H = (  0, -hz, hy)
  // (L =) M = E × n = (  0,  ez,-ey)
  if(IN_RT)
  {
    for(int j=sub_xbm; j<=sub_xtp; j++)
      for(int k=sub_xbk; k<=sub_xft; k++){
        CALC_COEFF(rt, j+subInfo_s.OFFSET_Y, k+subInfo_s.OFFSET_Z);
        int w = field_subIndex(sub_rt, j, k);
        SUB_EH_IN_YZ(w, ey,ez, hy,hz);
        Ny -= hz*coef;      Nz += hy*coef;
        Ly += ez*coef;      Lz -= ey*coef;
      }
  }
  
  //左の面 n=(-1, 0, 0)
  // (N =) J = n × H = (  0, hz,-hy)
  // (L =) M = E × n = (  0,-ez, ey)
  if(IN_LT)
  {
    for(int j=sub_xbm; j<=sub_xtp; j++)
      for(int k=sub_xbk; k<=sub_xft; k++){
        CALC_COEFF(lt, j+subInfo_s.OFFSET_Y, k+subInfo_s.OFFSET_Z);
        int w = field_subIndex(sub_lt, j, k);
        SUB_EH_IN_YZ(w, ey, ez, hy,hz);
        Ny += hz*coef;      Nz -= hy*coef;
        Ly -= ez*coef;      Lz += ey*coef;
      }
  }
  
  double sx = cos(theta_rad)*cos(phi_rad);
  double sy = cos(theta_rad)*sin(phi_rad);
  double sz = -sin(theta_rad);
  double px = -sin(phi_rad);
  double py = cos(phi_rad);

  dcomplex Nth = sx*Nx + sy*Ny + sz*Nz;
  dcomplex Nph = px*Nx + py*Ny;
  dcomplex Lth = sx*Lx + sy*Ly + sz*Lz;
  dcomplex Lph = px*Lx + py*Ly;

  //Coefficianは呼び出し元でかける. 
  *Eth = Coeffician*(Z0*Nth - Lph); //宇野先生の本では Nth, Nphに-がついていた.
  *Eph = Coeffician*(Z0*Nph + Lth);
}

//戻り値はRankが0か否か
static bool sumToRank0(dcomplex *p, int size)
{
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  
  if( subInfo_s.Rank == 0 )
  {
    MPI_Status status;
    dcomplex *tmp = newDComplex(size);
    for(int i=1; i<subInfo_s.Nproc; i++)
    {
      MPI_Recv(tmp, size, MPI_C_DOUBLE_COMPLEX, i, 0, MPI_COMM_WORLD, &status);
      for(int j=0; j<size; j++)
        p[j] += tmp[j];
    }
    free(tmp);
    return true;
  }
  else {
    MPI_Send(p, size, MPI_C_DOUBLE_COMPLEX, 0, 0, MPI_COMM_WORLD);
    return false;
  }
}

void ntff3D_SubFrequency( dcomplex *Ex, dcomplex *Ey,dcomplex *Ez,
                       dcomplex *Hx, dcomplex *Hy, dcomplex *Hz)
{
  dcomplex Eth[181][360], Eph[181][360];
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();

  double k_s = field_getK_S();
  dcomplex Coeffician = I * k_s / (4*M_PI*R) * cexp(-I*k_s*R);

  double ToRad = M_PI/180.0;
  for(int theta=0 ;theta<=180; theta++)
  {
    for(int phi=0 ; phi<360 ; phi++)
    {
      subFrequencyNTFF(Ex, Ey, Ez, Hx, Hy, Hz, &Eth[theta][phi], &Eph[theta][phi], theta*ToRad, phi*ToRad, Coeffician);
    }
  }  

  FieldInfo fInfo = field_getFieldInfo();
  if(sumToRank0(Eth, 181*360))
  {
    char buf[512];
    sprintf(buf, "%dnm_Eth_str.txt",fInfo.lambda_nm);
    FILE *fpEth = openFile(buf);
    for(int theta=0; theta<=180;theta++)
      {
	for(int phi=0; phi<360; phi++)
	  {
	    fprintf(fpEth, "%.20lf ", cnorm(Eth[theta][phi]));
	  }    
	fprintf(fpEth, "\n");
      }
    fclose(fpEth);
  }
  
  if(sumToRank0(Eph, 181*360))
  {
    char buf[512];
    sprintf(buf, "%dnm_Eph_str.txt",fInfo.lambda_nm);
    FILE *fpEph = openFile(buf);
    for(int theta=0; theta<=180; theta++)
      {
	for(int phi=0; phi<360; phi++)
	  {
	    fprintf(fpEph, "%.20lf ", cnorm(Eph[theta][phi]));
	  }
	fprintf(fpEph, "\n");
      }
    fclose(fpEph);
  }  
}
