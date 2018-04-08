classdef Aniso < handle

    properties
        address
    end

    methods
        function this = Aniso(ds, qr, ks, as, sr, fn, fm)
            this.address = AnisoWrapper('new', ds, qr, ks, as, sr, fn, fm);
            cprintf('*Green', 'Initializing Aniso at %d.\n', this.address);
        end

        function delete(this)
            cprintf('*Red', '  Destroying Aniso at %d.\n', this.address);
            AnisoWrapper('delete', this.address);
        end
        
        function setCoeff(this, sigma_s, sigma_t)
            AnisoWrapper('setCoeff', this.address, sigma_s, sigma_t);
        end
        
        function node = getNodes(this)
            node = AnisoWrapper('getNodes', this.address);
        end
        
        function cache(this)
            AnisoWrapper('cache', this.address);
        end
        
        function ret = mapping(this, charge)
            ret = AnisoWrapper('mapping', this.address, charge);
        end
       
    end
end