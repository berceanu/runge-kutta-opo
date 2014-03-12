      Program filter_en
        USE nag_fft, ONLY: nag_fft_2d, nag_fft_trig
        USE global
		use subroutines
		implicit none

        call read 

        allocate(trig_x(2*Nx), trig_y(2*Ny))
        allocate(y_tot_0(Nx,Ny,Nt))
        allocate(y_enfilt(Nx,Ny))

        call nag_fft_trig(trig_x)
        call nag_fft_trig(trig_y)

        call import_evolution

		!file for checking momentum and energy conservation
        open(unit=25, file="en-peaks.dat", status='replace')
        
        !filter the 3 peaks in energy (integrating)
		call filter_peak(omega_p, omega_cut_p, 'pump')		
        call write_kx_max('pump')
        call write_peak_mom('pump')
        y_enfilt(:,:)=nag_fft_2d(y_enfilt(:,:),inverse=.true.,trig_m=trig_x,trig_n=trig_y)    
        call write_peak_spc('pump')

		call filter_peak(omega_s, omega_cut_s, 'signal')		
        call write_kx_max('signal')        
        call write_peak_mom('signal')
        y_enfilt(:,:)=nag_fft_2d(y_enfilt(:,:),inverse=.true.,trig_m=trig_x,trig_n=trig_y)    
        call write_peak_spc('signal')

		call filter_peak(omega_id, omega_cut_id, 'idler')		
        call write_kx_max('idler')        
        call write_peak_mom('idler')
        y_enfilt(:,:)=nag_fft_2d(y_enfilt(:,:),inverse=.true.,trig_m=trig_x,trig_n=trig_y)    
        call write_peak_spc('idler')
		
        close(25)

        deallocate(trig_x, trig_y)
        deallocate(y_tot_0)
        deallocate(y_enfilt)

      end Program filter_en
