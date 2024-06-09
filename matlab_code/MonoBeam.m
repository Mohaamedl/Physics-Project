
classdef MonoBeam < handle
    properties
        f   % frequency of beam (THz)
        Delta_x  % FWHM of the cross-sectional intensity (mm)
        I0  % beam intensity at its center (TW/m^2)
        eta_x  % ratio of sampling region to beam FWHM
        Nx  % number of points to sample
        
        dx  % sampling interval
        xs  % position array
        lambda   % wavelength
        Rp  % radius of curvature
        z   % current propagation distance
        z0  % initial propagation distance
        w0  % beam waist radius
        E   % electric field array
    end
    
    methods
        function obj = MonoBeam(f, Delta_x, I0, eta_x, Nx)
            obj.f = f;  % frequency of beam (THz)
            obj.Delta_x = Delta_x;  % FWHM of the cross-sectional intensity (mm)
            obj.I0 = I0;  % beam intensity at its center (TW/m^2)
            obj.eta_x = eta_x;  % ratio of sampling region to beam FWHM
            obj.Nx = Nx;  % number of points to sample
            
            % Calculating other properties
            obj.dx = obj.Delta_x / obj.Nx;
            obj.xs = circshift((-obj.Delta_x / 2:obj.dx:obj.Delta_x / 2 - obj.dx)', obj.Nx / 2);
            obj.lambda = 3e8 / obj.f;  % Speed of light in mm/s
            obj.Rp = inf;
            obj.z = 0;
            obj.z0 = 0;
            obj.w0 = obj.Delta_x / sqrt(2 * log(2));
            E0 = sqrt(2 * pi * obj.I0);  % Amplitude of electric field
            obj.E = E0 * exp(-(obj.xs .^ 2) / (obj.w0 ^ 2));
        end
        
        function rotate(obj, alpha)
            alpha = deg2rad(alpha);
            if obj.dx >= (obj.lambda / (2 * sin(abs(alpha)))) && alpha ~= 0
                error('Phase aliasing occurred upon rotation');
            end
            obj.addPhase(tan(alpha) * obj.xs);
        end
        
        function mask(obj, M)
            M0 = M(0);
            phi0 = atan2(imag(M0), real(M0));
            obj.E = obj.E .* M(obj.xs) .* exp(-1i * phi0);
        end
        
        function propagate(obj, Delta_z)
            Delta_z = Delta_z * 3e8;  % Speed of light in mm/s
            zR = pi * obj.w0 ^ 2 / obj.lambda;
            dw = obj.z - obj.z0;
            propagator = [isinf(obj.Rp), abs(dw + Delta_z) < 10 * zR];
            
            if propagator == 1
                obj.propagateP2P(Delta_z);
            elseif propagator == 2
                obj.propagateP2P(-dw);
                obj.propagateW2S(Delta_z + dw);
                obj.Rp = Delta_z + dw;
            elseif propagator == 3
                obj.propagateS2W(-dw);
                obj.propagateP2P(Delta_z + dw);
                obj.Rp = inf;
            else
                obj.propagateS2W(-dw);
                obj.propagateW2S(Delta_z + dw);
                obj.Rp = Delta_z + dw;
            end
        end
        
        function lens(obj, f)
            f = f * 3e8;  % Speed of light in mm/s
            dw = obj.nudge(obj.z - obj.z0, 0);
            zR = pi * obj.w0 ^ 2 / obj.lambda;
            w = obj.w0 * sqrt(1 + (dw / zR) ^ 2);
            
            if dw ~= 0
                R = obj.nudge(dw * (1 + (zR / dw) ^ 2), f);
                if R ~= f
                    R = 1 / (1 / R - 1 / f);
                else
                    R = inf;
                end
            else
                R = -f;
            end
            
            if isfinite(R)
                eta = obj.lambda * R / (pi * w ^ 2);
                dw = R / (1 + eta ^ 2);
                w0 = w / sqrt(1 + 1 / eta ^ 2);
            else
                dw = 0;
                w0 = w;
            end
            
            zR = pi * w0 ^ 2 / obj.lambda;
            propagator = [sinf(obj.Rp), abs(dw) < 10 * zR];
            Rp = inf * (propagator(2)) + dw * (~propagator(2));
            
            if propagator == 1
                a = 1 / f;
            elseif propagator == 2
                a = 1 / f + 1 / Rp;
            elseif propagator == 3
                a = 1 / f - 1 / obj.Rp;
            else
                a = 1 / f - 1 / obj.Rp + 1 / Rp;
            end
            
            obj.addPhase(-abs(obj.xs) .^ 2 * a / 2);
            obj.z0 = obj.z - dw;
            obj.w0 = w0;
            obj.Rp = Rp;
        end
        
        function phi_u = getPhase(obj, filter_)
            E = obj.getField();
            phi = atan2(imag(E), real(E));
            phi_u = unwrap(phi);
            phi_u = phi_u + phi(obj.Nx / 2 + 1) - phi_u(obj.Nx / 2 + 1);
            if filter_
                phi_u = phi_u .* (abs(E) > max(abs(E(:))) * 1e-6);
            end
        end
        
        function A = getAmplitude(obj)
            A = abs(obj.E);
            A = circshift(A, obj.Nx / 2);
        end
        
        function E = getField(obj)
            E = obj.E;
            E = circshift(E, obj.Nx / 2);
        end
        
        function I = getIntensity(obj)
            I = abs(obj.E) .^ 2 / (2 * pi);
            I = circshift(I, obj.Nx / 2);
        end
        
        function plotPhase(obj, title, filter_)
            phi_u = obj.getPhase(filter_);
            plot(circshift(obj.xs, obj.Nx / 2), phi_u);
            title(title);
            xlabel('Position (m)');
            ylabel('Phase');
            if isinf(obj.Rp)
                warning('The beam cross-section is planar.');
            else
                warning('The beam cross-section is spherical with radius of curvature %.2f meters.', obj.Rp / 3e8);
            end
            grid on;
        end
        
        function plotAmplitude(obj, title)
            A = obj.getAmplitude();
            plot(circshift(obj.xs, obj.Nx / 2), A);
            title(title);
            xlabel('Position (m)');
            ylabel('Amplitude');
            grid on;
        end
        
        function plotIntensity(obj, title)
            I = obj.getIntensity();
            xs = circshift(obj.xs, obj.Nx / 2);
            plot(xs, I);
            title(title);
            xlabel('Position (m)');
            ylabel('Intensity');
            grid on;
        end

    end
    
    methods (Access = private)
        function obj = addPhase(obj, Delta_z)
            obj.E = obj.E .* exp(-2 * pi * 1i * Delta_z / obj.lambda);
        end
        
        function obj = propagateP2P(obj, Delta_z)
            if abs(Delta_z) < 1e-6
                return;
            end
            
            obj.computeFftw('FFTW_FORWARD');
            inverse_squared = circshift(((0:obj.Nx - 1) - obj.Nx / 2) / (obj.Nx * obj.dx) .^ 2, obj.Nx / 2);
            obj.E = obj.E .* exp(1i * pi * obj.lambda * Delta_z * inverse_squared);
            obj.computeFftw('FFTW_BACKWARD');
            obj.z = obj.z + Delta_z;
        end
        
        function obj = propagateW2S(obj, Delta_z)
            obj.addPhase(obj.xs .^ 2 / (2 * Delta_z));
            if Delta_z >= 0
                obj.computeFftw('FFTW_FORWARD');
            else
                obj.computeFftw('FFTW_BACKWARD');
            end
            
            dx = obj.lambda * abs(Delta_z) / (obj.dx * obj.Nx);
            obj.xs = obj.xs .* (dx / obj.dx);
            obj.Delta_x = obj.Delta_x * (dx / obj.dx);
            obj.dx = dx;
            obj.z = obj.z + Delta_z;
        end
        
        function obj = propagateS2W(obj, Delta_z)
            dx = obj.lambda * abs(Delta_z) / (obj.dx * obj.Nx);
            obj.xs = obj.xs .* (dx / obj.dx);
            obj.Delta_x = obj.Delta_x * (dx / obj.dx);
            obj.dx = dx;
            if Delta_z >= 0
                obj.computeFftw('FFTW_FORWARD');
            else
                obj.computeFftw('FFTW_BACKWARD');
            end
            
            obj.addPhase(obj.xs .^ 2 / (2 * Delta_z));
            obj.z = obj.z + Delta_z;
        end
        
        function obj = computeFftw(obj, direction)
            obj.E = fft(obj.E, obj.Nx);
            if strcmp(direction, 'FFTW_FORWARD')
                obj.E = obj.E / sqrt(obj.Nx);
            else
                obj.E = obj.E * sqrt(obj.Nx);
            end
        end
        
        function val = nudge(~, x, y)
            if isclose(x, y, 'RelTol', 1e-6)
                val = y;
            else
                val = x;
            end
        end
    end
end
