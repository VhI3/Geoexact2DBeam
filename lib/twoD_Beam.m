% /*
% MIT License
% 
% Copyright (c) 2024 VHI3
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy of
% this software and associated documentation files (the "Software"), to deal in
% the Software without restriction, including without limitation the rights to use,
% copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
% Software, and to permit persons to whom the Software is furnished to do so,
% subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
% HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
% ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
% THE USE OR OTHER DEALINGS IN THE SOFTWARE.
% */
%%
% Class definition for a 2D beam structure.
classdef twoD_Beam
    % Define properties (variables) of the class.
    properties
        D               % Constitutive matrix
        Structure       % Type of structure (e.g., 'Lee Frame', 'William Toggle', 'Cantilever Beam')
        nx              % Number of elements in the x-direction
        ny              % Number of elements in the y-direction
        p               % Polynomial degree
        ne              % Total number of elements
        nodCor          % Coordinates of the nodes
        eleNod          % Elements nodes
        GDofs           % Global degrees of freedom
        nNod            % Total number of nodes
        eleAng          % Angles of elements
        eleLen          % Lengths of elements
        fixNods         % Fixed nodes
        presDofs        % Prescribed degrees of freedom
        actvDofs        % Active degrees of freedom
        nActvDofs       % Number of active degrees of freedom
        maxDeflection   % Maximum deflection
        pointOfLoad     % Point of load application
        deformedShape   % Coordinates of deformed shape
        diaplacement    % Displacement vector
        actvDofsDisplacement % Displacement of active degrees of freedom
        lamda_deflection    % Lambda deflection
    end

    % Define methods (functions) of the class.

    methods
        % Constructor method for the twoD_Beam class.
        function obj = twoD_Beam(Structure, nx, ny, p)
            % Initialize properties with the provided arguments.
            obj.p = p;                % Polynomial degree
            ne = nx + ny;             % Total number of elements (sum of elements in x and y directions)
            obj.ne = ne;              % Store total number of elements
            obj.nx = nx;              % Number of elements in x-direction
            obj.ny = ny;              % Number of elements in y-direction
            nodCor = zeros((ne*p)+1,2); % Initialize node coordinates matrix

            % Conditional structure initialization based on the type of structure.
            if strcmp(Structure,'Lee Frame')
                % Specific initialization for 'Lee Frame'.
                Lx = 300; Ly = 300;    % Define dimensions for Lee Frame
                % Set y-coordinates for nodes in Lee Frame
                nodCor(1:ny*p + 1,2) = linspace(0, Ly, ny*p + 1)';
                % Set x-coordinates for nodes in Lee Frame
                nodCor(ny*p + 1:ne*p + 1,1) = linspace(0, Lx, nx*p + 1)';
                nodCor(ny*p + 1:ne*p + 1,2) = Ly;
                obj.nodCor = nodCor;    % Store node coordinates in object
                nNod = size(nodCor,1);  % Calculate the total number of nodes
                obj.nNod = nNod;        % Store number of nodes in object
                % Define element angles and lengths for Lee Frame
                alphaX = 0; alphaY = pi/2;
                [obj.eleNod, obj.eleAng, obj.eleLen] = obj.get_eleNod_local(Lx, Ly, nx, ny, p, alphaX, alphaY);
                GDofs = 3 * nNod;       % Calculate global degrees of freedom
                obj.GDofs = GDofs;      % Store global degrees of freedom
                obj.pointOfLoad = 3 * p * (ny + 1) + 2; % Define point of load application

            elseif strcmp(Structure,'William Toggle')
                % Specific initialization for 'William Toggle'.
                Lx = 65.715; Ly = 0.98; % Define dimensions for William Toggle
                % Initialize and set node coordinates for William Toggle
                nodCor = zeros((ne*p)+1,2);
                nodCorX = linspace(0, Lx, ne*p + 1)';
                nodCorY1 = linspace(0, Ly, ny*p + 1)';
                nodCorY2 = linspace(Ly, 0, ny*p + 1)';
                nodCorY2 = nodCorY2(2:end);
                nodCor(1:ne*p + 1, 1) = nodCorX;
                nodCor(1:nx*p + 1, 2) = nodCorY1;
                nodCor(nx*p + 2:ne*p + 1, 2) = nodCorY2;
                obj.nodCor = nodCor;    % Store node coordinates
                nNod = size(nodCor, 1); % Calculate total number of nodes
                obj.nNod = nNod;        % Store number of nodes
                % Define element nodes, angles, and lengths for William Toggle
                [obj.eleNod, obj.eleAng, obj.eleLen] = obj.get_eleNod_local_wt(Lx, Ly, nx, ny, p);
                GDofs = 3 * nNod;       % Calculate global degrees of freedom
                obj.GDofs = GDofs;      % Store global degrees of freedom
                obj.pointOfLoad = 3 * (floor((ne + 1) / 2) + 1) - 1; % Define point of load application

            elseif strcmp(Structure,'Canteliver Beam')
                % Specific initialization for 'Cantilever Beam'.
                Lx = 300; Ly = 0;       % Define dimensions for Cantilever Beam
                % Set x-coordinates for nodes in Cantilever Beam
                nodCor(1:ny*p + 1, 1) = linspace(0, Ly, ny*p + 1)';
                nodCor(ny*p + 1:ne*p + 1, 1) = linspace(Ly, Lx + Ly, nx*p + 1)';
                nNod = size(nodCor, 1); % Calculate total number of nodes
                obj.nodCor = nodCor;    % Store node coordinates
                % Define element angles and lengths for Cantilever Beam
                alphaX = 0; alphaY = 0;
                [obj.eleNod, obj.eleAng, obj.eleLen] = obj.get_eleNod_local(Lx, Ly, nx, ny, p, alphaX, alphaY);
                GDofs = 3 * nNod;       % Calculate global degrees of freedom
                obj.GDofs = GDofs;      % Store global degrees of freedom
                obj.pointOfLoad = 3 * obj.eleNod(end) - 1; % Define point of load application
                obj.nNod = nNod;        % Store number of nodes

            else
                % Error handling for unknown structure type.
                error('The Type of structure is unknown based on the given name!');
            end

        end


        % Method to get the constitutive matrix of the beam.
        function D = get_Constitutive_Matrix_GeoExact(obj)
            % Return the constitutive matrix stored in the object's properties.
            D = obj.D;
        end

        % Method to get the nodal coordinates.
        function nodCor = get_nodCor(obj)
            % Return the nodal coordinates matrix stored in the object's properties.
            nodCor = obj.nodCor;
        end

        % Method to set the constitutive matrix for the beam.
        function obj = set_Constitutive_Matrix_GeoExact(obj, D)
            % Ensure the input matrix D is a 3x3 matrix. If not, throw an error.
            if size(D, 1) ~= 3 || size(D, 2) ~= 3
                error('Matrix D must be a 3x3 matrix.');
            end
            % Assign the input matrix D to the constitutive matrix property of the object.
            obj.D = D;
        end


        % Method to get the element nodes, angles, and lengths for a local beam structure.
        function [eleNod,eleAng,eleLen] = get_eleNod_local(~, Lx, Ly, nx, ny, p, alphaX, alphaY)
            % Calculate the total number of elements.
            ne = nx + ny;
            % Initialize counters and loop through each element to assign nodes, angles, and lengths.
            counter = 1;
            % Loop for elements along the y-axis.
            for i = 1:ny
                for j = 1:p+1
                    eleNod(i,j) = counter; % Assign element node
                    counter = counter + 1;
                end
                counter = counter - 1; % Adjust counter at the end of each element row
                eleAng(i,1) = alphaY;  % Assign angle for y-axis elements
                eleLen(i,1) = Ly/ny;   % Assign length for y-axis elements
            end
            % Loop for elements along the x-axis.
            for i = (ny+1):ne
                for j = 1:p+1
                    eleNod(i,j) = counter; % Assign element node
                    counter = counter + 1;
                end
                counter = counter - 1; % Adjust counter at the end of each element row
                eleAng(i,1) = alphaX;  % Assign angle for x-axis elements
                eleLen(i,1) = Lx/nx;   % Assign length for x-axis elements
            end
        end


        % Method to get element nodes, angles, and lengths for 'William Toggle' structure.
        function [eleNod,eleAng,eleLen] = get_eleNod_local_wt(~, Lx, Ly, nx, ny, p)
            % Calculate effective lengths for elements based on the structure's geometry.
            Lx1 = sqrt((Lx/2)^2 + Ly^2); % Effective length for one part of the structure
            Lx2 = sqrt((Lx/2)^2 + Ly^2); % Effective length for the other part of the structure
            ne = nx + ny;                % Total number of elements
            counter = 1;
            % Assign nodes, angles, and lengths for elements along the y-axis.
            for i = 1:ny
                for j = 1:p+1
                    eleNod(i,j) = counter; % Assign element node
                    counter = counter + 1;
                end
                counter = counter - 1;
                eleAng(i,1) = atan(2*Ly/Lx); % Calculate and assign angle for the element
                eleLen(i,1) = Lx1/ny;        % Assign length for the element
            end
            % Assign nodes, angles, and lengths for elements along the x-axis.
            for i = (ny+1):ne
                for j = 1:p+1
                    eleNod(i,j) = counter; % Assign element node
                    counter = counter + 1;
                end
                counter = counter - 1;
                eleAng(i,1) = -atan(2*Ly/Lx); % Calculate and assign angle for the element
                eleLen(i,1) = Lx2/nx;         % Assign length for the element
            end
        end

        % Method to get element nodes from the object properties.
        function eleNod = get_eleNode(obj)
            % Return the element nodes stored in the object's properties.
            eleNod = obj.eleNod;
        end

        function obj = set_dirichLetBC(obj, option)
            eleNod = obj.eleNod;
            if option == 1
                % Both sides are simply supported. It means, the displacement on
                % both directions, x and y, are zero but the rotation is free.
                % Suitable for Lee Frame!
                obj.fixNods = [eleNod(1) eleNod(end)];
                obj.presDofs =[eleNod(1) eleNod(1)+1 3*eleNod(end)-2 3*eleNod(end)-1];

                %
            elseif option == 2
                % Both sides are fixed. It means, the displacement on both
                % directions, x and y, are zero and the rotation is also zero.
                % Suitable for William Toggle!
                obj.fixNods = [eleNod(1) eleNod(end)];
                obj.presDofs = [eleNod(1) eleNod(1)+1 eleNod(1)+2 3*eleNod(end)-2 3*eleNod(end)-1 3*eleNod(end)];
                %
            elseif option == 3
                % One side is fixed. It displements and rotations are zero
                % on one end.
                obj.fixNods = [eleNod(1)];
                obj.presDofs =[eleNod(1) eleNod(1)+1 eleNod(1)+2];
                %
            else
                error('The given option for boundary condition is not valid.');
            end
            obj.actvDofs = setdiff((1:obj.GDofs)',[obj.presDofs]);
            obj.nActvDofs = obj.GDofs - size(obj.presDofs,2);
        end



        % Method to solve the non-linear equations of the beam structure.
        function obj = solveNonLinearEQ(obj, beam_Theory, maxDeflection, nonLinearSolver)
            % Set the maximum deflection property of the object.
            obj.maxDeflection = maxDeflection;

            % Determine the appropriate functions for stiffness and internal force
            % based on the specified beam theory.
            if strcmp(beam_Theory,'Geometrically Exact')
                % For 'Geometrically Exact' beam theory, use specific functions.
                stiffnessFunc = sprintf("formReducedTangentStiffness%d", 1);
                internlForceFunc = sprintf("formReducedInternalForce%d", 1);

            elseif strcmp(beam_Theory,'Geometrically Moderate Rotation')
                % For 'Geometrically Moderate Rotation' beam theory, use specific functions.
                stiffnessFunc = sprintf("formReducedTangentStiffness%d", 2);
                internlForceFunc = sprintf("formReducedInternalForce%d", 2);

            elseif strcmp(beam_Theory,'Bernoulli')
                % For 'Bernoulli' beam theory, use specific functions.
                stiffnessFunc = sprintf("formReducedTangentStiffness%d", 3);
                internlForceFunc = sprintf("formReducedInternalForce%d", 3);

            elseif strcmp(beam_Theory,'2nd_Order')
                % For '2nd Order' beam theory, use specific functions.
                stiffnessFunc = sprintf("formReducedTangentStiffness%d", 4);
                internlForceFunc = sprintf("formReducedInternalForce%d", 4);

            else
                % If the beam theory is not recognized, throw an error.
                error('The beam theory is not valid.');
            end

            activeDof = obj.nActvDofs;     % Number of active degrees of freedom
            aDof = obj.actvDofs;           % Active degrees of freedom
            iq_C = zeros(obj.GDofs,1);     % Initialize load vector
            iq_C(obj.pointOfLoad) = -1;    % Assign load at the point of load application
            iq = iq_C(aDof);               % Load vector for active degrees of freedom
            Lpoint = find(iq<0);           % Find the index of the loaded point in active DOFs
            if strcmp(nonLinearSolver,'Newton-Raphson')

                %% Newton-Raphson Method Implementation
                tol = 1e-3;                 % Tolerance for convergence
                newton = 10000;             % Maximum number of Newton-Raphson iterations
                maxiter = 100;              % Maximum number of iterations for inner loop

                al = 0.0;                   % Initialize arc-length parameter
                a = zeros(activeDof,1);  % Initialize vector for active degrees of freedom displacements
                f = zeros(activeDof,1);  % Initialize vector for internal forces
                df = zeros(activeDof,activeDof); % Initialize tangent stiffness matrix
                dfinv = zeros(activeDof,activeDof); % Initialize inverse of tangent stiffness matrix
                store = [0 0];              % Initialize storage for displacements and arc-length parameter

                % Main loop for Newton-Raphson iterations
                for i = 1:newton
                    % Check if displacement at loaded point is less than or equal to maximum deflection
                    if a(Lpoint) <= obj.maxDeflection
                        break;  % Exit loop if maximum deflection is reached
                    end

                    % Initialize incremental displacement vector
                    da = zeros(activeDof,1);
                    dl = 0.5; % Incremental load factor

                    % Update displacement and load factor
                    a = a + da;
                    al = al + dl;

                    % Calculate internal forces
                    F_int = feval(internlForceFunc, obj.GDofs, obj.presDofs, obj.ne, obj.eleNod, obj.eleAng, obj.eleLen, obj.D, obj.p, a);

                    f = F_int - al*iq; % Compute the difference between internal forces and external load

                    % Check for convergence
                    fcheck = sqrt(f'*f);
                    if fcheck < tol
                        % Store results if convergence criteria is met
                        store = [store; a(Lpoint) al];
                    else
                        % Iterative process to achieve convergence
                        iters = 0;
                        while fcheck > tol
                            iters = iters + 1;  % Increment iteration count

                            % Compute the tangent stiffness matrix and its inverse
                            df = feval(stiffnessFunc, obj.GDofs, obj.presDofs, obj.ne, obj.eleNod, obj.eleAng, obj.eleLen, obj.D, obj.p, a);
                            dfinv = inv(df);

                            % Update incremental displacement
                            da = -(dfinv*f);
                            a = a + da;

                            % Recalculate internal forces
                            F_int = feval(internlForceFunc, obj.GDofs, obj.presDofs, obj.ne, obj.eleNod, obj.eleAng, obj.eleLen, obj.D, obj.p, a);
                            f = F_int - al*iq;

                            % Recheck for convergence
                            fcheck = sqrt(f'*f);
                            if iters > maxiter
                                % Display error message and stop if convergence is not achieved
                                disp('Convergence cannot be achieved within the maximum number of iterations.');
                                error('Program stops');
                            end
                        end
                        % Store results after convergence
                        store = [store; a(Lpoint) al];
                    end
                end
                % Store the final active degrees of freedom displacement and lambda deflection
                obj.actvDofsDisplacement = a;
                obj.lamda_deflection = store;

                % Update the overall displacement vector
                obj.diaplacement = zeros(obj.GDofs,1);
                obj.diaplacement(obj.actvDofs,1) = obj.actvDofsDisplacement;

                % Implementing the Riks Method (arc-length method) for solving non-linear equations.
            elseif strcmp(nonLinearSolver,'Arch-Length')
                %% Initialization of variables for the Riks Method
                tol = 1e-4;                    % Tolerance for convergence
                psi = 0.1;                       % Scaling factor for the arc-length constraint
                dll = 0.01;                    % Incremental arc-length
                U = zeros(obj.GDofs,1);        % Initialize displacement vector




                a = zeros(activeDof,1);        % Initialize vector for active displacements
                f = zeros(activeDof,1);        % Initialize vector for internal forces
                df = zeros(activeDof,activeDof); % Initialize tangent stiffness matrix
                dfinv = zeros(activeDof,activeDof); % Initialize inverse of tangent stiffness matrix
                dls = zeros(2,1);              % Initialize arc-length increments
                dao = zeros(activeDof,1);      % Initialize change in active displacements
                al = 0.0;                      % Initialize the arc-length parameter
                store = [0 0];                 % Initialize storage for displacements and arc-length parameter
                riks = 20000;                  % Maximum number of Riks iterations
                maxiter = 30;                  % Maximum number of iterations for convergence in each Riks step

                %% Riks Method main loop
                for i = 1:riks
                    % Check for maximum deflection criterion
                    if a(Lpoint) <= obj.maxDeflection
                        break;
                    end

                    % Initialize variables for each iteration
                    da = zeros(activeDof,1);
                    dab = zeros(activeDof,1);
                    dat = zeros(activeDof,1);
                    dda1 = zeros(activeDof,1);
                    dda2 = zeros(activeDof,1);
                    dda = zeros(activeDof,1);
                    dalfa = zeros(activeDof,1);
                    dl = 0;

                    % Calculate tangent stiffness matrix and its inverse
                    df = feval(stiffnessFunc, obj.GDofs, obj.presDofs, obj.ne, obj.eleNod, obj.eleAng, obj.eleLen, obj.D, obj.p, a+da);
                    dfinv = inv(df);
                    dat = dfinv*iq;

                    % Calculate arc-length increments
                    [ddl1,ddl2] = arc(da,dab,dat,dl,iq,psi,dll);
                    dda1 = dab + ddl1 * dat;
                    dda2 = dab + ddl2 * dat;

                    % Determine the correct arc-length increment
                    det0 = det(df);
                    if sign(det0) == sign(ddl1)
                        dda = dda1;
                        ddl = ddl1;
                    else
                        dda = dda2;
                        ddl = ddl2;
                    end
                    dalfa = da + dda;
                    dlamda = dl + ddl;

                    % Calculate internal forces
                    F_int = feval(internlForceFunc, obj.GDofs, obj.presDofs, obj.ne, obj.eleNod, obj.eleAng, obj.eleLen, obj.D, obj.p, a+dalfa);
                    f = F_int - (al+dlamda)*iq;

                    % Check for convergence
                    fcheck = sqrt(f'*f);
                    if fcheck < tol
                        % Update variables if convergence is achieved
                        a = a + dalfa;
                        al = al + dlamda;
                        store = [store; a(Lpoint) al];
                        dao = dalfa;
                        dlo = dlamda;
                    else
                        % Iterative process for convergence
                        iters = 0;
                        while fcheck > tol
                            % Increment iteration count
                            iters = iters + 1;
                            da = dalfa;
                            dl = dlamda;

                            % Recalculate internal forces
                            F_int = feval(internlForceFunc, obj.GDofs, obj.presDofs, obj.ne, obj.eleNod, obj.eleAng, obj.eleLen, obj.D, obj.p, a + da);
                            f = F_int - (al + dl)*iq;

                            % Recalculate tangent stiffness matrix and its inverse
                            df = feval(stiffnessFunc, obj.GDofs, obj.presDofs, obj.ne, obj.eleNod, obj.eleAng, obj.eleLen, obj.D, obj.p, a + da);
                            dfinv = inv(df);
                            dab = -(dfinv*f);
                            dat = (dfinv*iq);

                            % Recalculate arc-length increments
                            [ddl1, ddl2] = arc(da, dab, dat, dl, iq, psi, dll);
                            dda1 = dab + ddl1 * dat;
                            dda2 = dab + ddl2 * dat;
                            det0 = det(df);

                            % Determine the correct arc-length increment
                            daomag = (dao'*dao);
                            if daomag == 0.
                                if sign(dl + ddl1) == sign(det0)
                                    dda = dda1;
                                    ddl = ddl1;
                                else
                                    dda = dda2;
                                    ddl = ddl2;
                                end
                            else
                                aux1 = ((da + dda1)'*dao);
                                aux2 = ((da + dda2)'*dao);
                                aux3 = dlamda * (dl + ddl1) * (iq'*iq);
                                aux4 = dlamda * (dl + ddl2) * (iq'*iq);
                                dot1 = aux1 + psi^2* aux3;
                                dot2 = aux2 + psi^2* aux4;
                                if dot1 > dot2
                                    dda = dda1;
                                    ddl = ddl1;
                                else
                                    dda = dda2;
                                    ddl = ddl2;
                                end
                            end
                            if ddl1 == ddl2
                                dda = dda1;
                                ddl = ddl1;
                            end
                            dalfa = da + dda;
                            dlamda = dl + ddl;

                            % Recalculate internal forces
                            F_int = feval(internlForceFunc, obj.GDofs, obj.presDofs, obj.ne, obj.eleNod, obj.eleAng, obj.eleLen, obj.D, obj.p, a + dalfa);
                            f = F_int - (al + dlamda)*iq;
                            fcheck = norm(f);

                            % Check iteration limit
                            if iters > maxiter
                                iters = maxiter + 1;
                                break;
                            end
                        end
                        % Adjust arc-length increment based on iterations
                        dll = dll*1.5/iters;
                        if iters > maxiter
                            % Display error message and stop if convergence is not achieved
                            disp('Convergence cannot be achieved within the maximum number of iterations.');
                            error('Program stops');
                        else
                            % Update variables after convergence
                            a = a + dalfa;
                            al = al + dlamda;
                            store = [store; a(Lpoint) al];
                            dao = dalfa;
                            dlo = dlamda;
                        end
                    end
                end
                % Store the final active degrees of freedom displacement and lambda deflection
                obj.actvDofsDisplacement = a;
                obj.lamda_deflection = store;

                % Update the overall displacement vector
                obj.diaplacement = zeros(obj.GDofs,1);
                obj.diaplacement(obj.actvDofs,1) = obj.actvDofsDisplacement;
            else
                % Throw an error if the provided solver name is not valid.
                error('The method for non-linear solver is not valid.');
            end





        end

        function plot_deformedShape(obj, scale, titleStr,nn)
            U = obj.diaplacement(1:3:3*obj.nNod);
            W = obj.diaplacement(2:3:3*obj.nNod);
            obj.deformedShape = [obj.nodCor(:,1) + scale*U obj.nodCor(:,2) + scale*W];
            % figure
            subplot(4,2,nn)
            % draw2D(obj.nodCor ,obj.eleNod,obj.p, 'r-.')
            p1 = plot(obj.nodCor(:,1),obj.nodCor(:,2),'r-.');
            % set(p1,'DisplayName','Non-deformed Shape')
            axis off
            hold on
            % draw2D(obj.deformedShape ,obj.eleNod,obj.p,'b-.')
            p2 = plot(obj.deformedShape(:,1),obj.deformedShape(:,2),'b-.');
            % set(p2,'DisplayName','Geometrically Exact Deformed Shape')
            axis square
            title("Lee Frame")
            % legend
            sgtitle(titleStr)
        end

        function plot_lamda_deflection(obj,nn)
            store = obj.lamda_deflection;
            subplot(4,2,nn)
            a1 = plot(store(:,1),store(:,2),'r-');
            M1 = 'Geometically Exact Beam Theory';
            set(gca,'Xdir','reverse')
            ylabel({'$\lambda$'}, 'interpreter','latex')

            xlabel({'Deflection'})
            title("Deflection Lambda")
            grid on
            grid minor
            axis square
            hold off
        end

    end

end
