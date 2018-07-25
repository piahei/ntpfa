classdef PressureOilWaterModelNTPFAopt < PressureOilWaterModel
%  subclass to build model for the NTPFA method using optimization in
%  the decomposition of the conormal vector. 
%
% SYNOPSIS: 
%   model = PressureOilWaterModelNTPFAopt(G,rock,fluid);
% 
% PARAMETERS: 
%   G    - Grid structure as described by grid_structure.
%
%   rock - Rock data structure with valid field `perm`.  The permeability
%          is assumed to be in measured in units of metres squared (m^2).
%          Use function `darcy` to convert from darcies to m^2, e.g::
%
%                 perm = convertFrom(perm, milli*darcy)
%
%          if the permeability is provided in units of millidarcies.
%
%          The field rock.perm may have ONE column for a scalar
%          permeability in each cell, TWO/THREE columns for a diagonal
%          permeability in each cell (in 2/3 D) and THREE/SIX columns for a
%          symmetric full tensor permeability.  In the latter case, each
%          cell gets the permeability tensor::
%
%                 K_i = [ k1  k2 ]      in two space dimensions
%                       [ k2  k3 ]
%
%                 K_i = [ k1  k2  k3 ]  in three space dimensions
%                       [ k2  k4  k5 ]
%                       [ k3  k5  k6 ]
%   fluid - Fluid ADI data structure given by initSimpleADIFluid()
%
% RETURNS: 
%   model - NTPFAopt model for further use in incompSinglePhaseNTPFA
%
% SEE ALSO:
%   'pressureOilWaterModelNTPFA','computeGeometry','makeRock',
%   'initSimpleADIFluid','incompSinglePhaseNTPFA'.

    properties
        fluxType
        fixPointIteration
    end
    
    methods
        function model = PressureOilWaterModelNTPFAopt(G, rock, fluid, varargin)
            if isempty(fluid)
                fluid = initSimpleADIFluid();
            end
            
            model = model@PressureOilWaterModel(G, rock, fluid);
            model.incTolPressure = 1e-8;
            model.fixPointIteration = false;
            model = merge_options(model, varargin{:});
            model.fluxType = 'ntpfa';
            model.operators.collocationSet = getCollactionSetOPT(G, model.rock);
            
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = pressureEquationOilWaterNTPFA(state0, state, model,...
                dt, ...
                drivingForces,...
                varargin{:});
            
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            [state, report] = updateState@PressureOilWaterModel(model, state, problem, dx, drivingForces);
        end
        
        function  [convergence, values, names] = checkConvergence(model, problem)
            [convergence, values, names] = checkConvergence@PressureOilWaterModel(model, problem);
        end
        
        function [kgrad, model] = getPermGradient(model, p, p0, forces, transMult)
            % Get gradient operator for K grad(p)
            switch lower(model.fluxType)
                case 'tpfa'
                    T = model.operators.T.*transMult;
                    C = model.operators.C;
                    Ct = bsxfun(@times, C, T);
                    kgrad = @(x) -Ct*x;
                case 'mpfa'
                    error('NotImplemented');
                case 'ntpfa'
                    assert(isfield(model.operators, 'collocationSet'));
                    if model.fixPointIteration
                        p = double(p);
                    end
                    N = model.operators.N;
                    T = computeNonLinearTransForOpt(model.G, model.operators.collocationSet,p,N);
                    kgrad = @(x) T{2}.*x(N(:,2)) - x(N(:,1)).*T{1};
                otherwise
                    error(['Unknown discretization ''', model.fluxType, '''']);
            end
        end
    end
    
end