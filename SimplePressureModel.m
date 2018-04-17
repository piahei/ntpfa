classdef SimplePressureModel < ReservoirModel
    properties
        % Increment tolerance for pressure. Computes convergence in
        % pressure as the reduction in increments (scaled by the min/max
        % pressure of the reservoir)
        incTolPressure
        % Boolean indicating if increment tolerance is being used
        useIncTol
    end
    
    methods
        function model = SimplePressureModel(G, rock, fluid, varargin)
            model = model@ReservoirModel(G, rock, fluid);
            % Reasonable defaults
            model.incTolPressure = 1e-3;
            model.useIncTol = true;
            model = merge_options(model, varargin{:});
            model.water = true;
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            assert(isempty(drivingForces.W), 'Wells not supported');
            
            op = model.operators;
            % Pressure at current timestep
            p = model.getProps(state, 'pressure');
            % Pressure at previous timestep (not currently used)
            p0 = model.getProps(state0, 'pressure');

            % Set up AD and tell solver that pressure is the variable
            p = initVariablesADI(p);
            primaryVars = {'pressure'};

            % Get transmissibility and gradient
            T = op.T;
            dp = op.Grad(p);
            
            % Compute flux
            v = -T.*dp;
            
            % Take divergence
            peq = op.Div(v);
            
            eqs = {peq};
            names = {'water'};
            types = {'cell'};
            
            % Get source terms/bc. Create dummy values for saturation,
            % mobility and rho since this is single-phase without gravity
            s = ones(model.G.cells.num, 1);
            sat = {s};
            mob = {s};
            rho = {s};
            
            [eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                             {p}, sat, mob, rho, ...
                                                                             {}, {}, ...
                                                                             drivingForces);
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            p0 = state.pressure;
            [state, report] = updateState@ReservoirModel(model, state, problem, dx, drivingForces);
            range = max(p0) - min(p0);
            if range == 0
                range = 1;
            end
            state.dpRel = (state.pressure - p0)./range;
        end
        
        function  [convergence, values, names] = checkConvergence(model, problem)
            [convergence, values, names] = checkConvergence@PhysicalModel(model, problem);
            if ~isnan(problem.iterationNo) && model.useIncTol
                if problem.iterationNo  > 1
                    values(1) = norm(problem.state.dpRel, inf);
                else
                    values(1) = inf;
                end
                convergence = [values(1) < model.incTolPressure, values(2:end) < model.nonlinearTolerance];
                names{1} = 'Delta P';
            end
        end
    end
end
