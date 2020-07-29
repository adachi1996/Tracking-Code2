module map_track
  use const_data
  implicit none
  contains
  subroutine input_map
    integer :: i,r_num,th_num,z_num,all_num                                   !各データ数

    !初期設定関係
    select case (select_FFAG) !VFFAかHFFAかで使用する静止質量エネルギーを変えている
      case(1)
        m_m0 = mp0
      case(2)
        m_m0 = m0
    end select
    h        = m_th_h*(pi/180.0) !刻み幅の設定値[deg.]から[rad.]への変換
    pth      = (m_T**2.0 + 2.0*m_T*m_m0)**0.5 !初期運動量の計算
    particle = (/m_t0 ,m_th0 ,m_r0+dr ,m_z0+dz ,0.0d0 ,pth ,0.0d0 ,0.0d0 ,0.0d0 ,0.0d0/) !初期入射条件の計算
    !particle = (/m_t0 ,m_th0 ,m_r0+dr ,m_z0+dz,&
    !             -0.000040932662490143934d0,0.20610407862073235d0,0.000023481998702962754d0,0.0d0 ,0.0d0 ,0.0d0/)

    !ファイルの読み込み
    open(18,file=trim(path_name)//trim(file_name1)//trim(file_name2)//'.table', status='old')           !ファイルを開く
    read (18,*) z_num ,r_num,th_num                                            !1行目の読み込み。各データ数を保存(足立Ver)
    !read (18,*) th_num,z_num,r_num                                            !1行目の読み込み。各データ数を保存(和賀Ver)
    print *, 'Import for ',trim(path_name)//trim(file_name1)//trim(file_name2)//'.table'

    all_num = r_num * th_num * z_num                                          !全データ数の計算
    print *, all_num,r_num,th_num,z_num

    do i = 2 , 8                                                              !いらない行数を飛ばすループ
      read (18,*)
    end do

    select case (select_func)
      case(1)
        call coodinate(all_num,r_num,th_num,z_num)   !読み込み＆座標変換用サブルーチンを呼ぶ
        call map_tracking(r_num , th_num , z_num)    !軌道計算コードの呼び出しサブルーチンを呼ぶ
      case(2)
        call coodinate(all_num,r_num,th_num,z_num)   !読み込み＆座標変換用サブルーチンを呼ぶ
        call closed_orbit_H(r_num , th_num , z_num)  !HFFAの閉軌道計算サブルーチンを呼ぶ
      case(3)
        call coodinate(all_num,r_num,th_num,z_num)   !読み込み＆座標変換用サブルーチンを呼ぶ
      case(4)
        call coodinate(all_num,r_num,th_num,z_num)   !読み込み＆座標変換用サブルーチンを呼ぶ
        call closed_orbit_V(r_num , th_num , z_num)  !VFFAの閉軌道計算その1サブルーチンを呼ぶ
      case(5)
        call coodinate(all_num,r_num,th_num,z_num)   !読み込み＆座標変換用サブルーチンを呼ぶ
        call closed_orbit_V2(r_num , th_num , z_num) !VFFAの閉軌道計算その2サブルーチンを呼ぶ
      case(6)
        call coodinate(all_num,r_num,th_num,z_num)   !読み込み＆座標変換用サブルーチンを呼ぶ
        call closed_orbit_V3(r_num , th_num , z_num) !VFFAの閉軌道計算その3サブルーチンを呼ぶ
      case(7)
        call coodinate(all_num,r_num,th_num,z_num)   !読み込み＆座標変換用サブルーチンを呼ぶ
        call closed_orbit_V4(r_num , th_num , z_num) !VFFAの複数エネルギー軌道計算用サブルーチンを呼ぶ
    end select
    close(18)
  end subroutine

  !=====================================================[サブルーチン①:磁場マップ読み込み]=====================================================
  subroutine coodinate(all_num , r_num , th_num , z_num)
  double precision :: dat_x   , dat_y    , dat_z   , dat_Bx , dat_By , dat_Bz !読み込んだデータ格納用
  double precision :: dat_r   , dat_th   ,           dat_Br , dat_Bth         !座標変換したデータ格納用
  double precision :: temp_r  , temp_th  , temp_z                             !一時保管用
  integer          :: all_num , r_num    , th_num  , z_num                    !各データ数
  integer          :: count_r , count_th , count_z , i  ,count                !カウント用変数

  !読み込んだデータを格納する配列のポインター宣言
  double precision, pointer, dimension(:,:) :: temp_B_data
  double precision, pointer, dimension(:)   :: temp_r_data, temp_th_data, temp_z_data
  allocate( temp_r_data(r_num) , temp_th_data(th_num) , temp_z_data(z_num) , temp_B_data(all_num,3))

  count_r = 1     ; count_z = 1     ; count_th = 1     ; count = 1            !カウント用変数の初期化
  temp_r  = 1.0d5 ; temp_z  = 1.0d5 ; temp_th  = 1.0d5                        !一時保管用変数の初期化

  if (select_func == 3) then
    !open(15,file='r_th_z_20190504_200250300_original.csv', status='replace')                                                  !ファイル作成
    open(30,file='r_' //trim(file_name1)//trim(file_name2)//'.csv', status='replace')                                                  !ファイル作成
    open(31,file='th_'//trim(file_name1)//trim(file_name2)//'.csv', status='replace')                                                  !ファイル作成
    open(32,file='z_' //trim(file_name1)//trim(file_name2)//'.csv', status='replace')                                                  !ファイル作成
    open(33,file='B_' //trim(file_name1)//trim(file_name2)//'.csv', status='replace')                                                  !ファイル作成
    print *, 'make file!!'
  end if
  do i = 1 , all_num                                                          !全データを読み込むためのループ
    read (18,*) dat_x , dat_y , dat_z , dat_Bx , dat_By , dat_Bz              !読み込み
    dat_r   = 1.0d-4*(nint(1.0d2*(dat_x**2.0d0 + dat_y**2.0d0)**0.5))         !r[m]への座標変換
    dat_z   = 1.0d-2*dat_z                                                    !z[m]への単位変換
    dat_th  = 1.0d-4*(nint(1.0d4*(m_th0+ATAN(dat_y/dat_x)*180.0d0/pi)))         !角度の計算[deg.]
    dat_Br  = 1.0d-4*( dat_Bx*COS(dat_th*pi/180.0d0) + dat_By*SIN(dat_th*pi/180.0d0))              ![gauss]から単位変換[T]
    dat_Bth = 1.0d-4*(-dat_Bx*SIN(dat_th*pi/180.0d0) + dat_By*COS(dat_th*pi/180.0d0))              ![gauss]から単位変換[T]
    dat_Bz  = 1.0d-4*dat_Bz                                                                         ![gauss]から単位変換[T]

    if (count_r  <=  r_num .and. temp_r  /=  dat_r ) then                     !rのデータ数だけ格納
      temp_r_data(count_r)   = dat_r     ;     count_r  = count_r  + 1     ;     temp_r = dat_r
      write(30,*) dat_r
    end if

    if (count_z  <=  z_num .and. temp_z  /=  dat_z ) then                     !zのデータ数だけ格納
      temp_z_data(count_z)   = dat_z     ;     count_z  = count_z  + 1     ;     temp_z = dat_z
      write(32,*) dat_z
    end if

    if (count_th <= th_num .and. temp_th /= dat_th ) then                     !thのデータ数だけ格納
      temp_th_data(count_th) = dat_th    ;     count_th = count_th + 1     ;     temp_th = dat_th
      write(31,*) dat_th
    end if

    if (i == count*all_num/10) then
      print *, 'Now Loading ...',count*10,'%'
      count = count + 1
    end if

    temp_B_data(i,1) = dat_Br  !Br 格納
    temp_B_data(i,2) = dat_Bth !Bth格納
    temp_B_data(i,3) = dat_Bz  !Bz 格納
    if (select_func == 3) then
      !write(15,*) dat_r ,',', dat_th ,',', dat_z ,',', dat_Br ,',', dat_Bth ,',', dat_Bz
      write(33,*) dat_Br ,',', dat_Bth ,',', dat_Bz
    end if
  end do

  if (select_func == 3) then
    !close(15)
    close(30)
    close(31)
    close(32)
    close(33)
  end if

  r_data  => temp_r_data            !r のポインタ設定
  th_data => temp_th_data           !thのポインタ設定
  z_data  => temp_z_data            !z のポインタ設定
  B_data  => temp_B_data            !磁場のポインタ設定

  r_min  = r_data(1)                !rの最小値
  th_min = th_data(1)               !thの最小値
  z_min  = z_data(1)                !zの最小値
  drt_r  = r_data(2)  -  r_data(1)  !Δr
  drt_th = th_data(2) - th_data(1)  !Δth
  drt_z  = z_data(2)  -  z_data(1)  !Δz
  end subroutine

  !=====================================================[サブルーチン②:ルンゲクッタ呼び出し＆保存]=====================================================
  subroutine map_tracking(r_num , th_num , z_num)
    double precision :: temp_dth , temp_th , temp_deg , d_pth , d_T
    integer          :: r_num , th_num , z_num
    integer          :: count_cell , count_dth
    integer(8)       :: count_th
    open(17,file=save_name, status='replace')                                                  !ファイル作成
    write (17,*) 't[s]',',','th[deg]',',','r[m]',',','z[m]',',',' &
                 Pr[MeV/c]',',','Pth[MeV/c]',',','Pz[MeV/c]',',','Br[T]',',','Bth[T]',',','Bz[T]'   !保存データの名前書き込み

    temp_dth  = nint(m_dth/m_th_h)                                                                      !保存する角度の整数化
    temp_th   = 0                                                                                   !初期化
    count_cell= nint(sym_th/m_th_h)
    count_dth = nint(m_dth/m_th_h)
    count_th  = 0

    d_T   = (exp(3.0*0.001) - 1.0)
    d_pth = (particle(5)**2.0 + particle(6)**2.0 + particle(7)**2.0)**0.5!(m_T**2.0 + 2.0*m_T*m_m0)**0.5

    do while (particle(2) < max_deg)                                                                !任意の周回数回るまで計算するループ
      temp_deg = nint(particle(2)/m_th_h)                                                             !角度の整数化
      !print *, particle(2),',',particle(3),',',particle(4)
      !if (nint(mod(temp_deg , temp_dth)) == 0) then                                                 !任意の角度の整数倍の時にtrue
      if (mod(count_th,count_dth) == 0) then
        write (17,'(e16.8,a,e16.8,a,e16.8,a,e16.8,a,e16.8,a,e16.8,a,e16.8,a,e16.8,a,e16.8,a,e16.8)')&
                     particle(1),',',particle(2),',',particle(3),',', &
                     particle(4),',',particle(5),',',particle(6),',', &
                     particle(7),',',particle(8),',',particle(9),',',particle(10)                   !データを逐次保存
        !write (*,*) particle(2),particle(3),particle(4)
        !if (temp_th /= int(particle(2)/360)) then
        if (mod(count_th,nint(360/m_th_h)) == 0) then
          !d_pth = (particle(5)**2.0 + particle(6)**2.0 + particle(7)**2.0)**0.5
          !particle(6) = particle(6) + 0.0007148428!d_pth*d_T
          print *, nint(particle(2)/360),'Turns'!,'T=', d_pth,  d_pth*d_T                                                        !1周ごとに表示
          temp_th = int(particle(2)/360)                                                           !更新
        end if
      end if
      !theta = mod(temp_deg,sym_th/m_th_h)*m_th_h                                                          !1Cell内での角度の更新
      theta = (mod(count_th,count_cell))*m_th_h
      call map_cal(r_num , th_num , z_num)                                                          !RKを呼ぶための外部サブルーチンを呼ぶ
      particle(2) = particle(2) + m_th_h                                                              !角度の更新
      count_th = count_th + 1

      if (abs(m_z0 - particle(4)) > cut_z .or. check == 1) then                              !zの振幅が設定した範囲の外に出たら計算を終了する
        check = 0
        exit
      end if
    end do
    close(17)                                                                                       !ファイルを閉じる
  end subroutine

  !=====================================================[サブルーチン②:ルンゲクッタ計算読み出しと値の更新]=====================================================
  subroutine map_cal(r_num , th_num , z_num)
      double precision :: kl(6) , ks(6)  , k1(6) , k2(6) , k3(6) , k4(6) , k(6)                     !ルンゲクッタ用変数入れ
      double precision :: P
      integer          :: r_num , th_num , z_num

      kl = (/particle(1),particle(3),particle(4),particle(5),particle(6),particle(7)/)              !粒子の情報の借り入れ子
      ks = kl                                                                                       !粒子の情報の借り入れ子2

      P      = (particle(5)**2.0+particle(6)**2.0+particle(7)**2.0)**0.5                            !運動量(Ps的な感じ)
      r_beta =    P/(m_m0**2.0+P**2.0)**0.5                                                          !非速度
      gamma  = 1.0/(1.0-r_beta**2.0)**0.5                                                           !ローレンツ因子

      call map_B(r_num , th_num , z_num) !磁場読み込み
      particle = (/ks(1),particle(2),ks(2),ks(3),ks(4),ks(5),ks(6),particle(8),particle(9),particle(10)/)
      !call map_B(r_num , th_num , z_num) !磁場読み込み
      k1 = map_RK(ks)                    !一回目
      ks = kl + 0.5*k1                   !更新

      particle = (/ks(1),particle(2),ks(2),ks(3),ks(4),ks(5),ks(6),particle(8),particle(9),particle(10)/)
      !call map_B(r_num , th_num , z_num) !磁場読み込み
      k2 = map_RK(ks)                    !二回目
      ks = kl + 0.5*k2                   !更新

      particle = (/ks(1),particle(2),ks(2),ks(3),ks(4),ks(5),ks(6),particle(8),particle(9),particle(10)/)
      !call map_B(r_num , th_num , z_num) !磁場読み込み
      k3 = map_RK(ks)                    !三回目
      ks = kl + 1.0*k3                   !更新

      particle = (/ks(1),particle(2),ks(2),ks(3),ks(4),ks(5),ks(6),particle(8),particle(9),particle(10)/)
      !call map_B(r_num , th_num , z_num) !磁場読み込み
      k4 = map_RK(ks)                    !四回目

      k = (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
      kl = kl + k                        !最終更新

      particle = (/kl(1),particle(2),kl(2),kl(3),kl(4),kl(5),kl(6),particle(8),particle(9),particle(10)/)
  end subroutine

  !=====================================================[サブルーチン③:磁場補完部]=====================================================
  subroutine map_B(r_num , th_num , z_num)
    double precision :: sector_V , sum_V
    double precision :: temp_z
    integer :: i , j , k , n_r , n_th  , n_z
    integer :: r_num , th_num  , z_num , data_num

    i = int((particle(3)     - r_min )/drt_r ) + 1       !配列番号は1から始まるので + 1
    j = int((theta           - th_min)/drt_th) + 1       !配列番号は1から始まるので + 1
    k = int(abs((particle(4) - z_min )/drt_z)) + 1       !配列番号は1から始まるので + 1


    if (i >= r_num .or. j >= th_num .or. k >= z_num .or. i <= 0 .or. k <= 0) then !範囲外に出たときプログラムを止める
      !print *,i,j,k                                    !計算で求めた配列番号の表示
      !print *,r_num , th_num  , z_num                  !配列番号の最大値表示
      temp_th_max = particle(2)
      check = 1
      return
    end if

    !変数の初期化
    sum_V = 0.0d0
    particle( 8) =  0.0d0   ;   particle( 9) =  0.0d0   ;   particle(10) =  0.0d0

    do n_r  = 0 , 1    !r の変数
      do n_th = 0 ,1   !thの変数
        do n_z  = 0 ,1 !z の変数
          !体積の計算関係
          temp_z   = abs(abs(particle(4))-abs(z_data(k+n_z)))
          sector_V = 0.5*abs(particle(3)**2.0-r_data(i+n_r)**2.0) * abs(theta - th_data(j+n_th))*pi/180.0 * temp_z
          sum_V    = sum_V + sector_V
          !対応する配列番号の計算(体積と対角の座標を求める) 例 => (i + (1-n_r) -1)
          data_num = (i-n_r)*z_num*th_num+(k-n_z)*th_num+(j+1-n_th)
          !磁場を加算
          particle( 8) = particle( 8) + B_data( data_num ,1)*sector_V
          particle( 9) = particle( 9) + B_data( data_num ,2)*sector_V
          particle(10) = particle(10) + B_data( data_num ,3)*sector_V
        end do
      end do
    end do

    !磁場の更新
    if (particle(4) >= 0) then
      particle( 8) =  particle( 8)/sum_V   ;   particle( 9) =  particle( 9)/sum_V   ;   particle(10) =  particle(10)/sum_V
    else
      particle( 8) = -particle( 8)/sum_V   ;   particle( 9) = -particle( 9)/sum_V   ;   particle(10) =  particle(10)/sum_V
    end if
  end subroutine

  !=====================================================[サブルーチン④:ルンゲクッタ計算部]=====================================================
  function map_RK(ks) result(klf)                                                          !微分方程式の計算
    double precision :: ks(6),klf(6)
    klf(1) = h*(ks(2)*gamma*m_m0/ks(5))/c                                                   !t
    klf(2) = h*(ks(2)*ks(4)/ks(5))                                                         !r
    klf(4) = h*(q*(1.0d-6*c*ks(2)/ks(5))*(ks(5)*particle(10) - ks(6)*particle( 9)) + ks(5))  !pr
    klf(5) = h*(q*(1.0d-6*c*ks(2)/ks(5))*(ks(6)*particle( 8) - ks(4)*particle(10)) - ks(4))  !pth
    if (particle(4) == 0.0 .and. particle(7) == 0.0) then                                  !実数型のカスを消すためのif文
      klf(3) = 0.0d0
      klf(6) = 0.0d0
    else
      klf(3) = h*(ks(2)*ks(6)/ks(5))                                                       !z
      klf(6) = h*(q*(1.0d-6*c*ks(2)/ks(5))*(ks(4)*particle( 9) - ks(5)*particle( 8)))        !pz
    end if
  end function
    !=====================================================[サブルーチン⑤:複数エネルギーの閉軌道導出部(水平用)]=====================================================
  subroutine closed_orbit_H(r_num , th_num , z_num)
    double precision :: min_r , ini_r
    integer :: r_num , th_num , z_num , int_T ,i

    ini_r = m_r0
    do i = 20, 21
      T = dble(i)*0.001
      min_r  = 1.0d5                                                                    !初期化
      pth = (T**2.0 + 2.0*T*m_m0)**0.5
      particle = (/0.0d0, m_th0, ini_r, m_z0, 0.0d0, pth, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)           !粒子の初期情報を更新
      do while (abs(min_r) >= 1.0d-12)
        print *, min_r,ini_r
        ini_r  = particle(3)
        particle = (/0.0d0, m_th0, ini_r, m_z0, 0.0d0, pth, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)           !粒子の初期情報を更新
        do while (particle(2) < sym_th+m_th_h)                                                                !任意の周回数回るまで計算するループ
          theta = particle(2)                                                                        !1Cell内での角度の更新
          !print *, particle
          call map_cal(r_num , th_num , z_num)                                                          !RKを呼ぶための外部サブルーチンを呼ぶ
          particle(2) = particle(2) + m_th_h                                                              !角度の更新
        end do
        if (particle(2) >= sym_th) then
          min_r = (particle(3) - ini_r)
          if (min_r > 0) then
            ini_r = ini_r + min_r*0.5
          end if
        end if
      end do
      print *, T,'[MeV],r0=',ini_r , min_r
      int_T = int(T)

      write (filename1, '("closed_T=", i3.3, ".csv")') i ! ここでファイル名を生成している
     ! print *, trim(filename1)
      particle = (/0.0d0, m_th0, ini_r         , 0.0d0 , 0.0d0, pth, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)           !粒子の初期情報を更新
      call map_tracking2(r_num , th_num , z_num , 0)                                 !軌道計算コードの呼び出し

      write (filename2, '("1e-4_tune_T=", i3.3, ".csv")') i ! ここでファイル名を生成している
      !print *, trim(filename2)
      particle = (/0.0d0, m_th0, ini_r + 1.0d-3, 1.0d-3, 0.0d0, pth, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)           !粒子の初期情報を更新
      call map_tracking2(r_num , th_num , z_num , 1)                                 !軌道計算コードの呼び出し
    end do
  end subroutine

  !=====================================================[サブルーチン⑥:閉軌道の計算とチューン軌道計算]=====================================================
  subroutine map_tracking2(r_num , th_num , z_num , count)
    double precision :: temp_dth , temp_th , temp_deg
    double precision :: r , z
    integer          :: r_num , th_num , z_num , count
    if (count == 1) then
      open(20,file=trim(filename1), status='old')
      open(17,file=trim(filename2), status='replace')                                               !ファイル作成
      write (17,*) 't[s]',',','th[deg]',',','Δr[m]',',','Δz[m]',',','r[m]',',','z[m]'             !保存データの名前書き込み
    else
      open(17,file=trim(filename1), status='replace')                                               !ファイル作成
    end if

    temp_dth  = nint(m_dth/m_th_h)                                                                      !保存する角度の整数化
    temp_th   = 0                                                                                   !初期化

    do while (particle(2) < max_deg)                                                                !任意の周回数回るまで計算するループ
      !print *, particle(4)
      temp_deg = nint(particle(2)/m_th_h)                                                             !角度の整数化
      if (nint(mod(temp_deg , temp_dth)) == 0) then                                                 !任意の角度の整数倍の時にtrue
        if (count == 1) then
          read(20,*) r , z
          write (17,*) particle(1)   ,',', particle(2)    , ',' ,&
                       particle(3)-r ,',', particle(4) -z , ',' ,  r  ,',',  z                      !データを逐次保存
        else
          write (17,*) particle(3)   ,',', particle(4)
        end if
      end if
      theta = mod(temp_deg,sym_th/m_th_h)*m_th_h                                                          !1Cell内での角度の更新
      call map_cal(r_num , th_num , z_num)                                                          !RKを呼ぶための外部サブルーチンを呼ぶ
      particle(2) = particle(2) + m_th_h                                                              !角度の更新
    end do
    if (count == 1) then
      close(17)
      close(20,status='delete')
    else
      close(17)
    end if                                                                                 !ファイルを閉じる
  end subroutine
  !=====================================================[サブルーチン⑤:閉軌道導出(垂直用)]=====================================================
  subroutine closed_orbit_V(r_num , th_num , z_num)
    double precision :: temp_r1 , temp_r2 , ini_r                                              !ｒの情報の入れ物
    double precision :: temp_z1 , temp_z2 , ini_z                                              !ｚの情報の入れ物
    double precision :: temp_min                                                              !入り口と出口の差の入れ物
    double precision :: temp_deg                                                                    !角度用の入れ物
    double precision :: ini_dk , th_max
    integer          :: r_num , th_num , z_num
    integer          :: i,j,k,i_z                                                                                !整数変数

    ini_dk = dk
    do i_z = 0 , 1
      th_max = 0.0d0
      temp_r1 = m_r0
      temp_z1 = m_z0 + i_z*0.001
      temp_r2 = temp_r1
      temp_z2 = temp_z1
      temp_min = 10000.0d0
      write(*,*) "[r0,z0]=",temp_r1 , temp_z1
      do i = 1 , 7                                                                                   !小数点以下iまで計算するループ
        do j = -10 , 10                                                                               !ｒを変化させるためのループ
          ini_r = temp_r1 + dble(j)*dk                                                                !ｒの初期値を更新
          do  k = -10 , 10                                                                            !ｚを変化させるためのループ
            ini_z = temp_z1 + dble(k)*dk                                                           !ｚの初期値を更新
            particle = (/m_t0 ,m_th0, ini_r, ini_z, 0.0d0, pth, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)           !粒子の初期情報を更新
            do while (particle(2) <= 360.0d0)                                                         !1周回るまで計算するループ
              temp_deg = nint(particle(2)/m_th_h)                                                             !角度の整数化
              theta = mod(temp_deg,sym_th/m_th_h)*m_th_h                                                          !1Cell内での角度の更新
              call map_cal(r_num , th_num , z_num)                                                          !RKを呼ぶための外部サブルーチンを呼ぶ
              particle(2) = particle(2) + m_th_h                                                            !角度の更新
              if (abs(ini_z - particle(4)) > cut_z .or. check == 1) then                              !zの振幅が設定した範囲の外に出たら計算を終了する
                check = 0
                exit
              end if
            end do
            !print *, ini_r,ini_z,theta
            if (temp_min > ((ini_r-particle(3))**2. + (ini_z-particle(4))**2.)**0.5 &                 !より差が小さくなった時true
                                                  .and. nint(particle(2)) >= 360) then
              temp_min = ((ini_r-particle(3))**2. + (ini_z-particle(4))**2.)**0.5                     !差の更新
              temp_r2   = ini_r                                                                       !rの更新(仮)
              temp_z2   = ini_z                                                                       !zの更新(仮)
            end if
            if (temp_th_max > th_max .and. nint(temp_th_max) < 360) then
              th_max = temp_th_max
              temp_r_max = ini_r
              temp_z_max = ini_z
            end if
          end do !kのループ
        end do !jのループ
        dk      = dk*0.1                                                                              !刻み幅の更新
        temp_r1 = temp_r2                                                                             !rの更新(確定)
        temp_z1 = temp_z2                                                                             !zの更新(確定)
      end do !iのループ
      dk = ini_dk
      write(*,*) 'temp_min , th_max -->',temp_min , th_max
      write(*,*) 'r , z -->',temp_r2 , temp_z2
    end do !zのループ
  end subroutine

!=====================================================[サブルーチン⑤:角度も考慮した閉軌道導出(垂直用)]=====================================================
  subroutine closed_orbit_V2(r_num , th_num , z_num)
    double precision :: temp_r1  , temp_r2  , ini_r                                              !ｒの情報の入れ物
    double precision :: temp_z1  , temp_z2  , ini_z                                              !ｚの情報の入れ物
    double precision :: temp_pz1 , temp_pz2 , ini_pz                                             !Pｚの情報の入れ物
    double precision :: temp_pth1 , temp_pth2 , ini_pth                                             !Pｚの情報の入れ物
    double precision :: temp_min                                                              !入り口と出口の差の入れ物
    double precision :: temp_deg                                                                    !角度用の入れ物
    double precision :: ini_dk , th_max
    integer          :: r_num , th_num , z_num
    integer          :: i,j,k,l,i_z                                                                                !整数変数

    ini_dk = dk
    do i_z = 0 , 1
      th_max = 0.0d0
      temp_r1  = m_r0
      temp_z1  = m_z0 + i_z*0.001
      temp_pz1 = 0.0
      temp_r2  = temp_r1
      temp_z2  = temp_z1
      temp_pz2 = temp_pz1
      temp_min = 10000.0d0
      write(*,*) "[r0,z0]=",temp_r1 , temp_z1
      do i = 1 , 3                                                                                   !小数点以下iまで計算するループ
        do j = -10 , 10                                                                               !ｒを変化させるためのループ
          ini_r = temp_r1 + dble(j)*dk                                                                !ｒの初期値を更新
          do k = -10 , 10                                                                            !ｚを変化させるためのループ
            ini_z = temp_z1 + dble(k)*dk                                                           !ｚの初期値を更新
            do l = -10 , 10
              ini_pz   = temp_pz1 + dble(l)*dk
              ini_pth  = (pth**2.0 - ini_pz**2.0)**0.5

              particle = (/m_t0 ,m_th0, ini_r, ini_z, 0.0d0, ini_pth, ini_pz, 0.0d0, 0.0d0, 0.0d0/)           !粒子の初期情報を更新
              do while (particle(2) <= 360.0d0)                                                         !1周回るまで計算するループ
                temp_deg = nint(particle(2)/m_th_h)                                                             !角度の整数化
                theta = mod(temp_deg,sym_th/m_th_h)*m_th_h                                                          !1Cell内での角度の更新
                call map_cal(r_num , th_num , z_num)                                                          !RKを呼ぶための外部サブルーチンを呼ぶ
                particle(2) = particle(2) + m_th_h                                                            !角度の更新
                if (abs(ini_z - particle(4)) > cut_z .or. check == 1) then                              !zの振幅が設定した範囲の外に出たら計算を終了する
                  check = 0
                  exit
                end if
              end do
              !print *, ini_r,ini_z,theta
              if (temp_min > ((ini_r  -particle(3))**2. + (ini_z -particle(4))**2. + &
                              (ini_pth-particle(6))**2. + (ini_pz-particle(7))**2.)**0.5 &
                              .and. nint(particle(2)) >= 360) then                !より差が小さくなった時true
                 temp_min  = ((ini_r  -particle(3))**2. + (ini_z -particle(4))**2. + &
                              (ini_pth-particle(6))**2. + (ini_pz-particle(7))**2.)**0.5                     !差の更新
                 temp_r2   = ini_r                                                                       !rの更新(仮)
                 temp_z2   = ini_z                                                                       !zの更新(仮)
                 temp_pz2  = ini_pz                                                                       !zの更新(仮)
                 temp_pth2 = ini_pth
                 !write(*,*) th0,ini_r,ini_z
                 !write(*,*) particle(2),particle(3),particle(4)
                 !write(*,*) '==============================================',temp_min
               end if
               if (temp_th_max > th_max .and. nint(temp_th_max) < 360) then
                 th_max = temp_th_max
                 temp_r_max = ini_r
                 temp_z_max = ini_z
                 temp_pth_max = ini_pth
                 temp_pz_max  = ini_pz
                 write(*,*) ini_r , ini_z, ini_pth , ini_pz , temp_min
               end if
            end do !lのループ
          end do !kのループ
        end do !jのループ
        dk      = dk*0.1                                                                              !刻み幅の更新
        temp_r1 = temp_r2                                                                             !rの更新(確定)
        temp_z1 = temp_z2                                                                             !zの更新(確定)
        !write(*,*) i,temp_min,'dr,dz=',temp_r1-r0,temp_z1
      end do !iのループ
      dk = ini_dk
      write(*,*) 'temp_min , th_max -->',temp_min , th_max
      write(*,*) 'r   , z  -->',temp_r_max   , temp_z_max
      write(*,*) 'pth , pz -->',temp_pth2 , temp_pz2
    end do !zのループ
  end subroutine

  !=====================================================[サブルーチン⑤:角度も考慮した閉軌道導出(垂直用)]=====================================================
  subroutine closed_orbit_V3(r_num , th_num , z_num)
    double precision :: temp_r   , temp_r2    !rの情報の入れ物
    double precision :: temp_z   , temp_z2    !zの情報の入れ物
    double precision :: temp_pr  , temp_pr2   !prの情報の入れ物
    double precision :: temp_pz  , temp_pz2   !pzの情報の入れ物
    double precision :: temp_pth , temp_pth2  !pthの情報の入れ物
    double precision :: temp_deg                          !角度用の入れ物
    double precision :: temp_box                          !同じ結果を除くためのもの
    double precision :: max_r , max_z ,T_keV                         !同じ結果を除くためのもの
    integer          :: r_num , th_num , z_num            !各データ数用入れ子
    integer          :: i,j,k,set_num                        !ループ用整数変数
    integer          :: save_count,val_T ,i_r                          !保存名用整数変数
    integer          :: CO_z1,CO_z2,CO_z3
    integer          :: count_dth,count_th                !角度計算用整数変数
    character(100)   :: num_name                          !保存名用
    character(100)   :: add_name, temp_T                          !保存名用
    double precision, pointer, dimension(:) :: temp_sum   !平均値を求めるための配列
    double precision, pointer, dimension(:)   :: closed_z
    allocate( closed_z(20) )
    allocate( temp_sum(6) )                               !カウント数とr～pzまでの変数用

    set_num = 10
    temp_box   = 1.0d0
    save_count = 1
    count_dth  = nint(360.0/m_th_h) !保存する角度(360度で固定)
    do val_T = 1 , 25
      T_keV = dble(val_T)*2.0+20.0
      !print *, 'T =',nint(T_keV) ,' [keV]'
      pth = ((0.001*T_keV)**2.0 + 2.0*(0.001*T_keV)*m_m0)**0.5
      CO_z1 = 0
      CO_z2 = 0
      do CO_z3 = 1 , 20
        closed_z(CO_z3) = 0.0
      end do
    do i = 0 , 240      !zのループ
    do i_r = -6 , 6      !rのループ
      max_r = 0.0
      max_z = 0.0
      temp_pr2  = 0.0d0 !prの初期化
      temp_pz2  = 0.0d0 !pzの初期化
      temp_pth2 = pth   !pthの初期化
      temp_r2   = m_r0 + dble(i_r)*0.005  !rの初期化
      temp_z2   = m_z0 + dble(i  )*0.001   !zの初期化
      !if (mod(i,24) == 0) then
        !print *, '   x,y = ',temp_r2,temp_z2
      !end if
      do k = 0 , set_num            !入射条件の更新回数
        if (k >= set_num) then
          write (num_name, '("alpha=175_B0=130_m=3_", i3.3, "_")') save_count
          write (add_name, '("T=", i3.3, ".csv")') nint(T_keV)
          save_count = save_count + 1
          open(20,file=trim(num_name)//trim(add_name), status='replace')                                                  !ファイル作成
          write (20,*) 't[s]',',','th[deg]',',','r[m]',',','z[m]',',',' &
                        Pr[MeV/c]',',','Pth[MeV/c]',',','Pz[MeV/c]',',','Br[T]',',','Bth[T]',',','Bz[T]'   !保存データの名前書き込み
        end if
        temp_sum(1) = 0.0      !22.5度ごとのデータの総数
        temp_sum(2) = 0.0      !22.5度ごとのrの総和
        temp_sum(3) = 0.0      !22.5度ごとのzの総和
        temp_sum(4) = 0.0      !22.5度ごとのprの総和
        temp_sum(6) = 0.0      !22.5度ごとのpzの総和
        temp_sum(5) = 0.0      !22.5度ごとのpthの総和
        temp_r   = temp_r2     !初期値rの更新
        temp_z   = temp_z2     !初期値zの更新
        temp_pr  = temp_pr2    !初期値prの更新
        temp_pz  = temp_pz2    !初期値pzの更新
        temp_pth = temp_pth2   !初期値pthの更新
        count_th = 0           !カウント数の初期化
        !print *,k,max_r,max_z
        if (k <  set_num) particle = (/m_t0 ,m_th0, temp_r+max_r, temp_z+max_z, temp_pr, temp_pth, temp_pz, 0.0d0, 0.0d0, 0.0d0/) !粒子の初期情報を更新
        if (k == set_num) particle = (/m_t0 ,m_th0, temp_r      , temp_z      , temp_pr, temp_pth, temp_pz, 0.0d0, 0.0d0, 0.0d0/) !粒子の初期情報を更新
        max_r = 0.0
        max_z = 0.0
        do while (particle(2) <= (360.0d0*50.0))                                                    !1周回るまで計算するループ
          temp_deg = nint(particle(2)/m_th_h)                                                       !角度の整数化
          theta    = mod(temp_deg,sym_th/m_th_h)*m_th_h                                             !1Cell内での角度の更新
          if (mod(count_th,count_dth) == 0) then    !22.5度の倍数の時に入る
            temp_sum(1) = temp_sum(1) + 1           !データ数の計算(総和)
            temp_sum(2) = temp_sum(2) + particle(3) !rの総計
            temp_sum(3) = temp_sum(3) + particle(4) !zの総計
            temp_sum(4) = temp_sum(4) + particle(5) !Prの総計
            temp_sum(5) = temp_sum(5) + particle(6) !Pthの総計
            temp_sum(6) = temp_sum(6) + particle(7) !Pzの総計
            if (k == set_num) then
              write (20,'(e16.8,a,e16.8,a,e16.8,a,e16.8,a,e16.8,a,e16.8,a,e16.8,a,e16.8,a,e16.8,a,e16.8)')&
                           particle(1),',',particle(2),',',particle(3),',', &
                           particle(4),',',particle(5),',',particle(6),',', &
                           particle(7),',',particle(8),',',particle(9),',',particle(10)                   !データを逐次保存
            end if
            if (max_r < particle(3)) max_r = particle(3)
            if (max_z < particle(4)) max_z = particle(4)
          end if
          call map_cal(r_num , th_num , z_num)                                                      !RKを呼ぶための外部サブルーチンを呼ぶ
          particle(2) = particle(2) + m_th_h                                                        !角度の更新
          if (abs(temp_z - particle(4)) > cut_z .or. check == 1) then                               !zの振幅が設定した範囲の外に出たら計算を終了する
            check = 1 !kのループからも抜けるために値を更新
            exit      !軌道計算のループから抜ける
          end if
          count_th = count_th + 1 !角度用カウント数の更新
        end do
        if (check == 1) then      !任意のターン数(多分30周)回らなかった場合に入る
          check = 0               !計算終了条件用変数を初期状態に戻す
          if (k >= set_num) close(20)               !開いていたファイルを閉める
          exit                    !次のzに移るためにkのループから抜ける
        end if
        temp_r2   = temp_sum(2)/temp_sum(1) !rの更新(確定)
        temp_z2   = temp_sum(3)/temp_sum(1) !zの更新(確定)
        temp_pr2  = temp_sum(4)/temp_sum(1) !prの更新(確定)
        temp_pz2  = temp_sum(6)/temp_sum(1) !pzの更新(確定)
        temp_pth2 = (pth**2.0 - temp_pr2**2.0 - temp_pz2**2.0)**0.5
        !temp_pth2 = temp_sum(5)/temp_sum(1) !pthの更新(確定)
        max_r = (max_r - temp_r2)*0.1
        max_z = (max_z - temp_z2)*0.1
        if (k >= set_num) close(20)                                                                                !ファイルを閉じる


        if (CO_z1 >= 1) then
          do CO_z3 = 1, CO_z1
            !print *, '------->',CO_z3,closed_z(CO_z3)
            if (int(closed_z(CO_z3)*1d5) == int(temp_z2*1d5)) then
              CO_z2 = 1
            end if
          end do
        else if (CO_z1 == 0 .and. k == set_num) then
          CO_z1 = CO_z1 + 1
          closed_z(CO_z1) = temp_z2
        end if
        if (CO_z2 == 1) then
          !print *, 'k=',k
          exit
        end if
        if (k == set_num) then
          CO_z1 = CO_z1 + 1
          closed_z(CO_z1) = temp_z2
          !write(*,*) '------------------------------------------------------------'
          !write(*,*) temp_r2  , temp_z2
          !write(*,*) temp_pr2 , temp_pth2 , temp_pz
          !write(*,*) 'max_r_z=',max_r,max_z
          !write(*,*) 'diff=',((particle(3)-temp_r)**2.0 + (particle(4)-temp_z)**2.0)**0.5
          !write(*,*) '------------------------------------------------------------'
          write (temp_T,'(f0.1)') T_keV
          write(*,'(a,e16.8,e16.8,e16.8,e16.8,e16.8)') trim(temp_T),particle(3),particle(4),particle(5),particle(6),particle(7)
        end if
        CO_z2 = 0
      end do !kのループ
    end do !iのループ
    end do
    end do
  end subroutine

  !=====================================================[サブルーチン⑤:角度も考慮した閉軌道導出(垂直用)]=====================================================
  subroutine closed_orbit_V4(r_num , th_num , z_num)
  character(100)   :: read_file_name, temp_name
  integer          :: r_num   , th_num  , z_num     !各データ数用入れ子
  integer          :: i_num   , d_num               !各データ数用入れ子
  double precision :: temp_r  , temp_z  ,temp_T     !rの情報の入れ物
  double precision :: temp_Pr , temp_Pz , temp_Pth  !zの情報の入れ物

  open(40,file=trim(path_name2)//trim(file_name3)//'.txt', status='old')           !ファイルを開く
  read (40,*) read_file_name
  print *, 'Read file = '//read_file_name
  read (40,*) d_num
  do i_num = 1,d_num
    particle = (/m_t0 ,m_th0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
    read (40,*) temp_T,particle(3), particle(4), particle(5), particle(6), particle(7)
    m_z0 = particle(4)
    particle(3) = particle(3) + dr
    particle(4) = particle(4) + dz
    write (temp_name,'(f0.1)') temp_T
    save_name = trim(read_file_name)//'_T='//trim(temp_name)//'_dz=0.0005.csv'
    print *, nint(temp_T),save_name
    call map_tracking(r_num , th_num , z_num)
  end do
  close(40)
  end subroutine
end module map_track
