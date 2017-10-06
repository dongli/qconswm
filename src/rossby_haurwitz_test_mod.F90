module rossby_haurwitz_test_mod

  use mesh_mod
  use parallel_mod
  use params_mod
  use dycore_mod

  implicit none

  private

  public rossby_haurwitz_test_set_initial_condition

  integer, parameter :: R = 4
  real, parameter :: omg = 3.924e-6
  real, parameter :: gd0 = 8.0e3 * g

contains

  ! u = a ω (cosφ + R cosᴿ⁻¹φ sin²φ cosRλ - cosᴿ⁺¹φ sinφ cosRλ)
  !
  ! v = - a ω R cosᴿ⁻¹φ sinφ sinRλ
  !
  ! gd = gd0 + a² A(φ) + a² B(φ) cosRλ + a² C(φ) cos2Rλ
  !
  ! A(φ) = 1/2 ω (2 Ω + ω) cos²φ + 1/4 ω² cos²ᴿφ ((R + 1) cos²φ + (2 R² - R - 2) - 2 R² cos⁻²φ)
  ! B(φ) = 2 (Ω + ω) ω cosᴿφ ((R² + 2 R + 2) - (R + 1)² cos²φ) / (R + 1) / (R + 2)
  ! C(φ) = 1/4 ω² cos²ᴿφ ((R + 1) cos²φ - (R + 2))


  subroutine rossby_haurwitz_test_set_initial_condition()

    real a, b, c
    integer i, j

    write(6, *) '[Notice]: Use Rossby-Haurwitz wave initial condition.'

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        static%ghs(i,j) = 0.0
      end do
    end do

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        a = mesh%full_cos_lat(j)
        b = R * mesh%full_cos_lat(j)**(R - 1) * mesh%full_sin_lat(j)**2 * cos(R * mesh%half_lon(i))
        c = mesh%full_cos_lat(j)**(R + 1) * mesh%full_sin_lat(j) * cos(R * mesh%half_lon(i))
        state(1)%u(i,j) = radius * omg * (a + b - c)
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        a = R * mesh%half_cos_lat(j)**(R - 1) * mesh%half_sin_lat(j) * sin(R * mesh%full_lon(i))
        state(1)%v(i,j) = - radius * omg * a
      end do
    end do

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      a = 0.5 * omg * (2 * omega + omg) * mesh%full_cos_lat(j) + &
        0.25 * omg**2 * mesh%full_cos_lat(j)**(2 * R) * &
        ((R + 1) * mesh%full_cos_lat(j)**2 + 2 * R**2 - R - 2 - 2 * R**2 * mesh%full_cos_lat(j)**(-2))
      b = 2 * (omega + omg) * omg * mesh%full_cos_lat(j)**R * &
        (R**2 + 2 * R + 2 - (R + 1)**2 * mesh%full_cos_lat(j)**2) / (R + 1) / (R + 2)
      c = 0.25 * omg**2 * mesh%full_cos_lat(j)**(2 * R) * ((R + 1) * mesh%full_cos_lat(j)**2 - R - 2)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state(1)%gd(i,j) = gd0 + radius**2 * (a + b * cos(R * mesh%full_lon(i)) + c * cos(2 * R * mesh%full_lon(i)))
      end do
    end do

  end subroutine rossby_haurwitz_test_set_initial_condition

end module rossby_haurwitz_test_mod