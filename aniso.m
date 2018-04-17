classdef aniso < handle
    
    properties (Access = public)
        rte
        fem
        DSA
        f2rId
        r2fId
        intrId
        edgeId
        Diff
        Diff2
        sigma_s_node
        sigma_a_node
        sigma_t_node
        g
        N
        n
    end
    
    methods
        function obj = aniso(g, N)
            % generate RTE
            obj.rte = Aniso(128, 1, N, g, 10, 4, 20);
            nodes = obj.rte.getNodes();
            
            obj.n = size(nodes, 1);
            
            obj.N = N;
            obj.g = g;
            % n = 0 .. 2 * N - 2
            
            % generate mesh
            mesh = TriangleMesh();
            augmentedNodes = nodes';
            augmentedNodes = augmentedNodes(:);
            
            % insert the vertices at the boundary
            augmentedNodes = [[0 0 1 0 1 1 0 1]' ; augmentedNodes];
            mesh.set_points_tri(augmentedNodes);
            mesh.set_facets_tri([ 0 1 1 2 2 3 3 0]');
            mesh = mesh.build_tri();
            mesh = mesh.refine_tri('q34.0a0.5');
            mesh.getInfo_tri()
            
            % generate fem (edge is not necessary now)
            opt = struct('deg', 1, 'qdeg', 2, 'min_area', 4e-3,...
                'edge', [0 1 1 0; 0 0 1 1], 'mesh', mesh);
            obj.fem = femm(opt);
            
            % reorder all the nodes.
            obj.edgeId = unique(obj.fem.space.edges)';
            obj.intrId = setdiff(1:size(obj.fem.space.nodes, 2), obj.edgeId);
            fullNodes = obj.fem.space.nodes';
            intrNodes = fullNodes(obj.intrId, :);
            [~, intrOrd] = sortrows(intrNodes, [1 2]);
            obj.f2rId = obj.intrId(intrOrd);
            rIntrOrd(intrOrd) = 1:length(intrOrd);
            obj.r2fId = (rIntrOrd);            
            obj.Diff = [];
            obj.Diff2 = [];
        end
        
        % generate the finite element matrix for diffusion equation with
        % given function handles.
        function diffgen(obj, sigma_s, sigma_a)
            % sigma_s and sigma_a are given by function handles.
            % charge inputs are only evaluated on the interior nodes
            % augment with background value to the boundary.
            % charge should be augmented by 0 in general.
%             augCharge = zeros(size(obj.fem.space.nodes, 2), 1); 
%             augCharge(obj.f2rId) = charge;
            
            param_s = sigma_s(obj.fem.space.nodes)';
            param_a = sigma_a(obj.fem.space.nodes)';
            
            D = 0.5 ./(param_s + param_a + zeros(size(obj.fem.space.nodes, 2), 1));
            A = param_a + zeros(size(obj.fem.space.nodes, 2), 1);
            T = (param_s + param_a + zeros(size(obj.fem.space.nodes, 2), 1));
            
            qD = obj.mapping(D, obj.fem.space.elems, obj.fem.facet.ref');
            qA = obj.mapping(A, obj.fem.space.elems, obj.fem.facet.ref');
            qT = obj.mapping(T, obj.fem.space.elems, obj.fem.facet.ref');
            
            S = obj.fem.build('s', qD);
            M = obj.fem.build('m', qA);
            J = obj.fem.build('m', qT);
            E = obj.fem.build('e', 1, 'all');
            
            obj.Diff = (S + M + 0.5 * E);
            obj.Diff2 = (S + J + 0.5 * E);
            
        end
        
        function setCoeff(obj, sigma_s, sigma_a)
            nodes = obj.rte.getNodes();
            obj.sigma_s_node = sigma_s(nodes)' + zeros(size(nodes, 1), 1);
            obj.sigma_a_node = sigma_a(nodes)' + zeros(size(nodes, 1), 1);
            obj.sigma_t_node = obj.sigma_s_node + obj.sigma_a_node;
            
            obj.rte.setCoeff(obj.sigma_s_node, obj.sigma_t_node);
            
            tic;
            for i = 0 : (2 * obj.N - 2)
                obj.rte.cache(i);
            end
            toc;

        end
        
        function x = prec(obj, h)
            T = zeros(size(obj.fem.space.nodes, 2), 1);
            T(obj.f2rId) = h;
            qT = obj.mapping(T, obj.fem.space.elems, obj.fem.facet.ref');
            L = obj.fem.build('l', qT);
            z = obj.Diff \ (obj.Diff2 * L);
            x = zeros(size(h, 1), 1);
            x(obj.r2fId) = z(obj.intrId);
        end        
        
        function rhs = forward(obj, u)
            rhs = zeros(obj.N * obj.n, 1);
            
            for i = -(obj.N-1):0 
                iid = abs(i);
                
                for j = -(obj.N - 1):(obj.N - 1)
                    id = abs(j);
                   
                    theta = obj.rte.mapping(u( (id * obj.n + 1) : ((id+1) * obj.n)  ),  abs(i-j));
   
                    rhs(((iid) * obj.n+1 ): (iid+1) * obj.n) =...
                            rhs(((iid) * obj.n+1 ): (iid+1) * obj.n) + theta;
                   
                end
            end
        end
        
        function rhs = mforward(obj, u)
            rhs = zeros(obj.N * obj.n, 1);
            for i = -(obj.N - 1):0
                iid = abs(i);
               
                for j = -(obj.N - 1):(obj.N - 1)
                    id = abs(j);
                    
                    theta = obj.rte.mapping( obj.sigma_s_node .* u( (id * obj.n + 1) : ((id+1) * obj.n)  ),  abs(i-j));
                    
                    chi = (obj.g^id - obj.g^(obj.N)) / (1 - obj.g^(obj.N));
                    
                    rhs(((iid) * obj.n+1 ): (iid+1) * obj.n) =...
                            rhs(((iid) * obj.n+1 ): (iid+1) * obj.n) - chi * theta;
                   
                end
            end
        end
        
        
        function [y, flag, relres, iter, resvec] = solve(obj, charge, pr)
            rhs = obj.forward(charge);
          
            A = @(x)(x - obj.mforward(x));
            if pr == 1
                tic;
                [y, flag, relres, iter, resvec] = gmres(A, rhs, 400, 1e-15, 400, @obj.prec);
                toc;
            else
                tic;
                [y, flag, relres, iter, resvec] = gmres(A, rhs, 400, 1e-15, 400,[]);
                toc;
            end
                
        end
        
    end
    
    
    methods(Static)
        function [interpolate] = mapping(func, elems, trans_ref)
            numberofqnodes = size(trans_ref, 1);
            interpolate = zeros(numberofqnodes, size(elems, 2));
            for i = 1: size(elems, 2)
                interpolate(:, i) = trans_ref * func(elems(:, i));
            end
        end
        function r = normsq(v)
            r = sum(v.^2);
        end
    end
end

