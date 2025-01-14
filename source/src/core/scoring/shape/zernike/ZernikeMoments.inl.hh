namespace core {
    namespace scoring {
        namespace shape {
            namespace zernike {

                template<class VoxelT, class MomentT>
                inline typename ZernikeMoments<VoxelT, MomentT>::ComplexT
                ZernikeMoments<VoxelT, MomentT>::GetMoment(int _n, int _l, int _m) {
                    if (_m >= 0) {
                        return zernikeMoments_[_n][_l / 2][_m];
                    } else {
                        T sign;
                        if (_m % 2) {
                            sign = (T)(-1);
                        } else {
                            sign = (T) 1;
                        }
                        return sign * std::conj(zernikeMoments_[_n][_l / 2][abs(_m)]);
                    }
                }

            }
        }
    }
}