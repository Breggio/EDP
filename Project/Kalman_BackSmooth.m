function [Z_smoothed, SmoothErr_CovMatr] = Kalman_BackSmooth(Z_filtered,PredErr_CovMatr, FiltrErr_CovMatr, T)
    size_ = length( Z_filtered );
    TransMatr = [1 T 0 0; 0 1 0 0; 0 0 1 T; 0 0 0 1];
    SmoothErr_CovMatr=FiltrErr_CovMatr;
	Z_smoothed = Z_filtered;
    
    for i = size_-1 : -1 : 1
        A_coef=FiltrErr_CovMatr(:,:,i)*TransMatr.'/PredErr_CovMatr(:,:,i+1);
        SmoothErr_CovMatr(:,:,i) = FiltrErr_CovMatr(:,:,i) + ...
            A_coef*(SmoothErr_CovMatr(:,:,i+1) - PredErr_CovMatr(:,:,i))*A_coef.';
        Z_smoothed(:,i) = Z_filtered(:,i) + ...
            A_coef*(Z_smoothed(:,i+1)-TransMatr*Z_filtered(:,i));
    end
end
