# Original problems from sddpdams
# TODO:
# - factorize dynamics and costs with matrix from:
#       x_{t+1} = A x_t + B u_t + b_t
#       D x_t + E u_t <= g_t
#       cost: c_t u_t + d_t x_t

import JuMP

#------------------------------
# 1 dam valley
#------------------------------
if VALLEY=="dam1"
    const N_DAMS = 1
    # Damsvalley configuration:
    const STOCK_TARGET = [40 for i in 1:N_DAMS]
    const LOWER_DAM = [2]

    const X_BOUNDS = [(0, 80)];
    const U_BOUNDS = vcat([(0, 40)],
                          [(0, Inf) for i in 1:N_DAMS]);

    # Discretization of controls for assessment:
    DU = [10, 8]

    function dynamic(t, x, u, w)
        return [x[1] - u[1] + w[1] - u[1+N_DAMS]]
    end

    if LINEAR_COST
        function cost_t(t, x, u, w)
            return - COST[t] * sum(u[1:N_DAMS])
        end
    else
        function cost_t(t, x, u, w)
            return (- COST[t] * sum(u[1:N_DAMS]) +
                      EPSILON_U*sum(u[1:N_DAMS]'*u[1:N_DAMS]) +
                      EPSILON_Z*sum((u[1:N_DAMS]+u[N_DAMS+1:2*N_DAMS])'*
                                    (u[1:N_DAMS]+u[N_DAMS+1:2*N_DAMS])))
        end
    end

    function final_cost_dams(model, m)
        alpha = JuMP.getindex(m, :alpha)
        w = JuMP.getindex(m, :w)
        x = JuMP.getindex(m, :x)
        u = JuMP.getindex(m, :u)
        xf = JuMP.getindex(m, :xf)
        @JuMP.variable(m, z1[1:N_DAMS] >= 0)
        @JuMP.constraint(m, alpha == 0.)
        @JuMP.constraint(m, z1 .>= STOCK_TARGET - xf)
        @JuMP.objective(m, Min, model.costFunctions(model.stageNumber-1,x,u,w) +
                                alpha + COST_COEFF*sum(z1'*z1))
    end

#------------------------------
# 2 dams valley
#------------------------------
elseif VALLEY=="dam2"
    const N_DAMS = 2
    # Damsvalley configuration:
    const STOCK_TARGET = Int64[STOCK0 for i in 1:N_DAMS]
    const LOWER_DAM = [2, 3]
    const X_BOUNDS = [(STOCK_MIN, STOCK_MAX) for i in 1:N_DAMS];
    const U_BOUNDS = vcat([(CONTROL_MIN, CONTROL_MAX) for i in 1:N_DAMS], [(0., 200) for i in 1:N_DAMS]);
    # Discretization of controls for assessment:
    # Values for dam2  : DU = Int64[8, 8, 8, 8]
    # Values for dam2b : DU = Int64[40, 40, 40, 40]
    DU = Int64[8, 8, 8, 8]

    # Define dynamic of the dam:
    function dynamic(t, x, u, w)
        return [x[1] - u[1] + w[1] - u[3],
                x[2] - u[2] + u[1] + u[3] + w[2] - u[4]]
    end

    # Define cost corresponding to each timestep:
    if LINEAR_COST
        function cost_t(t, x, u, w)
        return - COST[t] * sum(u[1:N_DAMS])
        end
    else
        function cost_t(t, x, u, w)
            return (- COST[t] * sum(u[1:N_DAMS]) + EPSILON_U*sum(u[1:N_DAMS]'*u[1:N_DAMS]) +
                        EPSILON_Z*sum((u[1:N_DAMS]+u[N_DAMS+1:2*N_DAMS])'*(u[1:N_DAMS]+u[N_DAMS+1:2*N_DAMS])))
        end
    end

    # Build final cost as a LP model:
    # TODO: DRY this code
    function final_cost_dams(model, m)
        alpha = JuMP.getvariable(m, :alpha)
        w = JuMP.getvariable(m, :w)
        x = JuMP.getvariable(m, :x)
        u = JuMP.getvariable(m, :u)
        xf = JuMP.getvariable(m, :xf)
        @JuMP.variable(m, z1 >= 0)
        @JuMP.variable(m, z2 >= 0)
        @JuMP.constraint(m, alpha == 0.)
        @JuMP.constraint(m, z1 >= 40 - xf[1])
        @JuMP.constraint(m, z2 >= 40 - xf[2])
        @JuMP.objective(m, Min, model.costFunctions(model.stageNumber-1, x, u, w) + alpha + COST_COEFF*(z1*z1+z2*z2))
    end

#------------------------------
# 3 dams valley
#------------------------------
elseif VALLEY=="dam3"
    const N_DAMS = 3
    # Damsvalley configuration:
    const STOCK_TARGET = Int64[STOCK0 for i in 1:N_DAMS]
    const LOWER_DAM = [2, 3, 4]
    const X_BOUNDS = [(STOCK_MIN, STOCK_MAX) for i in 1:N_DAMS];
    const U_BOUNDS = vcat([(CONTROL_MIN, CONTROL_MAX) for i in 1:N_DAMS], [(0., 200) for i in 1:N_DAMS]);
    # Discretization of controls for assessment:
    # Values for dam3  : DU = Int64[8, 8, 8, 8, 8, 8]
    # Values for dam3b : DU = Int64[40, 40, 40, 40, 40, 40]
    DU = Int64[8, 8, 8, 8, 8, 8]

    # Define dynamic of the dam:
    function dynamic(t, x, u, w)
        return [x[1] - u[1] + w[1] - u[4],
                x[2] - u[2] + u[1] + u[4] + w[2] - u[5],
                x[3] - u[3] + u[2] + u[5] + w[3] - u[6]]
    end

    # Define cost corresponding to each timestep:
    if LINEAR_COST
        function cost_t(t, x, u, w)
        return - COST[t] * sum(u[1:N_DAMS])
        end
    else
        function cost_t(t, x, u, w)
            return (- COST[t] * sum(u[1:N_DAMS]) + EPSILON_U*sum(u[1:N_DAMS]'*u[1:N_DAMS]) +
                        EPSILON_Z*sum((u[1:N_DAMS]+u[N_DAMS+1:2*N_DAMS])'*(u[1:N_DAMS]+u[N_DAMS+1:2*N_DAMS])))
        end
    end

    # Build final cost as a LP model:
    # TODO: DRY this code
    function final_cost_dams(model, m)
        alpha = JuMP.getvariable(m, :alpha)
        w = JuMP.getvariable(m, :w)
        x = JuMP.getvariable(m, :x)
        u = JuMP.getvariable(m, :u)
        xf = JuMP.getvariable(m, :xf)
        @JuMP.variable(m, z1 >= 0)
        @JuMP.variable(m, z2 >= 0)
        @JuMP.variable(m, z3 >= 0)
        @JuMP.constraint(m, alpha == 0.)
        @JuMP.constraint(m, z1 >= 40 - xf[1])
        @JuMP.constraint(m, z2 >= 40 - xf[2])
        @JuMP.constraint(m, z3 >= 40 - xf[3])
        @JuMP.objective(m, Min, model.costFunctions(model.stageNumber-1, x, u, w) + alpha + COST_COEFF*(z1*z1+z2*z2+z3*z3))
    end

#------------------------------
# 4 dams valley
#------------------------------
elseif VALLEY=="dam4"
    const N_DAMS = 4
    # Damsvalley configuration:
    const STOCK_TARGET = [STOCK0 for i in 1:N_DAMS]
    const LOWER_DAM = [2, 3, 4, 5]
    const X_BOUNDS = [(STOCK_MIN, STOCK_MAX) for i in 1:N_DAMS];
    const U_BOUNDS = vcat([(CONTROL_MIN, CONTROL_MAX) for i in 1:N_DAMS], [(0., Inf) for i in 1:N_DAMS]);

    # Discretization of controls for assessment:
    DU = [8, 8, 8, 8, 8, 8, 8, 8]

    function dynamic(t, x, u, w)
        return [x[1] - u[1] + w[1] - u[5],
                x[2] - u[2] + u[1] + u[5] + w[2] - u[6],
                x[3] - u[3] + u[2] + u[6] + w[3] - u[7],
                x[4] - u[4] + u[3] + u[7] + w[4] - u[8]]
    end

    if LINEAR_COST
        function cost_t(t, x, u, w)
            return - COST[t] * (u[1] + u[2] + u[3] + u[4])
        end
    else
        function cost_t(t, x, u, w)
            return - COST[t] * (u[1] + u[2] + u[3] + u[4]) + EPSILON_U*(u[1]*u[1] + u[2]*u[2] + u[3]*u[3] + u[4]*u[4] ) +
                        EPSILON_Z*((u[1]+u[5])*(u[1]+u[5]) + (u[2]+u[6])*(u[2]+u[6]) + (u[3]+u[7])*(u[3]+u[7])+ (u[4]+u[8])*(u[4]+u[8]))
        end
    end

    function final_cost_dams(model, m)
        alpha = JuMP.getvariable(m, :alpha)
        w = JuMP.getvariable(m, :w)
        x = JuMP.getvariable(m, :x)
        u = JuMP.getvariable(m, :u)
        xf = JuMP.getvariable(m, :xf)
        @JuMP.variable(m, z1 >= 0)
        @JuMP.variable(m, z2 >= 0)
        @JuMP.variable(m, z3 >= 0)
        @JuMP.variable(m, z4 >= 0)
        @JuMP.constraint(m, alpha == 0.)
        @JuMP.constraint(m, z1 >= 40. - xf[1])
        @JuMP.constraint(m, z2 >= 40. - xf[2])
        @JuMP.constraint(m, z3 >= 40. - xf[3])
        @JuMP.constraint(m, z4 >= 40. - xf[4])
        @JuMP.objective(m, Min, model.costFunctions(model.stageNumber-1, x, u, w) + alpha + COST_COEFF*(z1*z1+z2*z2+z3*z3+z4*z4))
    end

#------------------------------
# 5 dams valley
#------------------------------
elseif VALLEY=="dam5"
    const N_DAMS = 5
    # Damsvalley configuration:
    const STOCK_TARGET = [STOCK0 for i in 1:N_DAMS]
    const LOWER_DAM = [3, 3, 4, 5, 6]
    const X_BOUNDS = [(STOCK_MIN, STOCK_MAX) for i in 1:N_DAMS];
    const U_BOUNDS = vcat([(CONTROL_MIN, CONTROL_MAX) for i in 1:N_DAMS], [(0., Inf) for i in 1:N_DAMS]);

    # Discretization of controls for assessment:
    DU = [8, 8, 8, 8, 8, 8, 8, 8, 8, 8]

    function dynamic(t, x, u, w)
        return [x[1] - u[1] + w[1] - u[6],
                x[2] - u[2] + w[2] - u[7],
                x[3] - u[3] + u[2] + u[7] + u[6] + u[1] + w[3] - u[8],
                x[4] - u[4] + u[3] + u[8] + w[4] - u[9],
                x[5] - u[5] + u[4] + u[9] + w[5] - u[10]]
    end

    if LINEAR_COST
        function cost_t(t, x, u, w)
            return - COST[t] * (u[1] + u[2] + u[3] + u[4] + u[5])
        end
    else
        function cost_t(t, x, u, w)
            return - COST[t] * (u[1] + u[2] + u[3] + u[4] + u[5]) + EPSILON_U*(u[1]*u[1] + u[2]*u[2] + u[3]*u[3] + u[4]*u[4] + u[5]*u[5]) +
                        EPSILON_Z*((u[1]+u[6])*(u[1]+u[6]) + (u[2]+u[7])*(u[2]+u[7]) + (u[3]+u[8])*(u[3]+u[8])+ (u[4]+u[9])*(u[4]+u[9])+ (u[5]+u[10])*(u[5]+u[10]))
        end
    end

    function final_cost_dams(model, m)
        alpha = JuMP.getvariable(m, :alpha)
        x = JuMP.getvariable(m, :x)
        u = JuMP.getvariable(m, :u)
        w = JuMP.getvariable(m, :w)
        xf = JuMP.getvariable(m, :xf)
        @JuMP.variable(m, z1 >= 0)
        @JuMP.variable(m, z2 >= 0)
        @JuMP.variable(m, z3 >= 0)
        @JuMP.variable(m, z4 >= 0)
        @JuMP.variable(m, z5 >= 0)
        @JuMP.constraint(m, alpha == 0.)
        @JuMP.constraint(m, z1 >= 40. - xf[1])
        @JuMP.constraint(m, z2 >= 40. - xf[2])
        @JuMP.constraint(m, z3 >= 40. - xf[3])
        @JuMP.constraint(m, z4 >= 40. - xf[4])
        @JuMP.constraint(m, z5 >= 40. - xf[5])
        @JuMP.objective(m, Min, model.costFunctions(model.stageNumber-1, x, u, w) + alpha + COST_COEFF*(z1*z1+z2*z2+z3*z3+z4*z4+z5*z5))
    end

#------------------------------
# 6 dams valley
#------------------------------
elseif VALLEY == "dam6"
    const N_DAMS = 6
    # Damsvalley configuration:
    const STOCK_TARGET = [40, 40, 40, 40, 40, 40]
    const LOWER_DAM = [2, 4, 4, 5, 6, 7]
    const X_BOUNDS = [(0, 80), (0, 80), (0, 80), (0, 80), (0, 80), (0, 80)];
    const U_BOUNDS = vcat([(0, 40), (0, 40), (0, 40), (0, 80), (0, 64), (0, 64)],
                            [(0, Inf) for i in 1:N_DAMS]);

    # Discretization of controls for assessment:
    DU = [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8]

    function dynamic(t, x, u, w)
        return [x[1] - u[1] + w[1] - u[7],
                x[2] - u[2] + u[1] + u[7] + w[2] - u[8],
                x[3] - u[3] + w[3] - u[9],
                x[4] - u[4] + u[3] + u[2] + u[8] + u[9] + w[4] - u[10],
                x[5] - u[5] + u[4] + u[10] + w[5] - u[11],
                x[6] - u[6] + u[5] + u[11] + w[6] - u[12]]
    end

    if LINEAR_COST
        function cost_t(t, x, u, w)
            return - COST[t] * (u[1] + u[2] + u[3] + u[4] + u[5] + u[6])
        end
    else
        function cost_t(t, x, u, w)
            return - COST[t] * (u[1] + u[2] + u[3] + u[4] + u[5] + u[6]) + EPSILON_U*(u[1]*u[1] + u[2]*u[2] + u[3]*u[3] + u[4]*u[4] + u[5]*u[5] + u[6]*u[6]) +
                        EPSILON_Z*((u[1]+u[6])*(u[1]+u[6]) + (u[2]+u[7])*(u[2]+u[7]) + (u[3]+u[8])*(u[3]+u[8])+
                        (u[4]+u[9])*(u[4]+u[9])+ (u[5]+u[10])*(u[5]+u[10])+ (u[6]+u[11])*(u[6]+u[11]))
        end
    end

    function final_cost_dams(model, m)
        alpha = JuMP.getvariable(m, :alpha)
        w = JuMP.getvariable(m, :w)
        x = JuMP.getvariable(m, :x)
        u = JuMP.getvariable(m, :u)
        xf = JuMP.getvariable(m, :xf)
        @JuMP.variable(m, z1 >= 0)
        @JuMP.variable(m, z2 >= 0)
        @JuMP.variable(m, z3 >= 0)
        @JuMP.variable(m, z4 >= 0)
        @JuMP.variable(m, z5 >= 0)
        @JuMP.variable(m, z6 >= 0)
        @JuMP.constraint(m, alpha == 0.)
        @JuMP.constraint(m, z1 >= STOCK_TARGET[1] - xf[1])
        @JuMP.constraint(m, z2 >= STOCK_TARGET[2] - xf[2])
        @JuMP.constraint(m, z3 >= STOCK_TARGET[3] - xf[3])
        @JuMP.constraint(m, z4 >= STOCK_TARGET[4] - xf[4])
        @JuMP.constraint(m, z5 >= STOCK_TARGET[5] - xf[5])
        @JuMP.constraint(m, z6 >= STOCK_TARGET[6] - xf[6])
        @JuMP.objective(m, Min, model.costFunctions(model.stageNumber-1, x, u, w) + alpha + COST_COEFF*(z1*z1+z2*z2+z3*z3+z4*z4+z5*z5+z6*z6))
    end

#------------------------------
# 8 dams valley
#------------------------------
elseif VALLEY == "dam8"
    const N_DAMS = 8
    # Damsvalley configuration:
    const STOCK_TARGET = [40 for i in 1:N_DAMS]
    const LOWER_DAM = [2, 6, 5, 5, 6, 7, 8, 9]
    const X_BOUNDS = [(0, 80),(0, 100),(0, 80),(0, 80),(0, 100),(0, 100),(0, 100),(0, 100)];
    const U_BOUNDS = vcat([(0, 40), (0, 40), (0, 40), (0, 40), (0, 80), (0, 120), (0, 120), (0, 120)],
                            [(0, Inf) for i in 1:N_DAMS]);

    # Discretization of controls for assessment:
    DU = [10, 10, 10, 10, 10, 20, 20, 20, 8, 8, 8, 8, 8, 8, 8, 8]

    function dynamic(t, x, u, w)
        return [x[1] - u[1] + w[1] - u[1+N_DAMS],
                x[2] - u[2] + u[1] + u[1+N_DAMS] + w[2] - u[2+N_DAMS],
                x[3] - u[3] + w[3] - u[3+N_DAMS],
                x[4] - u[4] + w[4] - u[4+N_DAMS],
                x[5] - u[5] + u[3] + u[4] + u[3+N_DAMS] + u[4+N_DAMS] + w[5] - u[5+N_DAMS],
                x[6] - u[6] + u[2] + u[5] + u[2+N_DAMS] + u[5+N_DAMS] + w[6] - u[6+N_DAMS],
                x[7] - u[7] + u[6] + u[6+N_DAMS] + w[7] - u[7+N_DAMS],
                x[8] - u[8] + u[7] + u[7+N_DAMS] + w[8] - u[8+N_DAMS]]
    end

    if LINEAR_COST
        function cost_t(t, x, u, w)
            return - COST[t] * sum(u[1:N_DAMS])
        end
    else
        function cost_t(t, x, u, w)
            return (- COST[t] * sum(u[1:N_DAMS]) + EPSILON_U*sum(u[1:N_DAMS]'*u[1:N_DAMS]) +
                        EPSILON_Z*sum((u[1:N_DAMS]+u[N_DAMS+1:2*N_DAMS])'*(u[1:N_DAMS]+u[N_DAMS+1:2*N_DAMS])))
        end
    end

    function final_cost_dams(model, m)
        alpha = JuMP.getvariable(m, :alpha)
        w = JuMP.getvariable(m, :w)
        x = JuMP.getvariable(m, :x)
        u = JuMP.getvariable(m, :u)
        xf = JuMP.getvariable(m, :xf)
        @JuMP.variable(m, z1[1:N_DAMS] >= 0)
        @JuMP.constraint(m, alpha == 0.)
        @JuMP.constraint(m, z1 .>= STOCK_TARGET - xf)
        @JuMP.objective(m, Min, model.costFunctions(model.stageNumber-1, x, u, w) +
                                alpha + COST_COEFF*sum(z1'*z1))
    end

#------------------------------
# 10 dams valley
#------------------------------
elseif VALLEY == "damX"
    const N_DAMS = 10
    # Damsvalley configuration:
    const STOCK_TARGET = [40 for i in 1:N_DAMS]
    const LOWER_DAM = [2, 6, 5, 5, 6, 8, 8, 9, 10, 11]
    const X_BOUNDS = [(0, 80),(0, 100),(0, 80),(0, 80),(0, 100),(0, 100),(0, 80),(0, 100),(0, 100),(0, 100)];
    const U_BOUNDS = vcat([(0, 40), (0, 40), (0, 40), (0, 40), (0, 80), (0, 120), (0, 40), (0, 140), (0, 140), (0, 140)],
                            [(0, Inf) for i in 1:N_DAMS]);

    # Discretization of controls for assessment:
    DU = [10, 10, 10, 10, 10, 20, 10, 20, 20, 20, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8]

    function dynamic(t, x, u, w)
        return [x[1] - u[1] + w[1] - u[1+N_DAMS],
                x[2] - u[2] + u[1] + u[1+N_DAMS] + w[2] - u[2+N_DAMS],
                x[3] - u[3] + w[3] - u[3+N_DAMS],
                x[4] - u[4] + w[4] - u[4+N_DAMS],
                x[5] - u[5] + u[3] + u[4] + u[3+N_DAMS] + u[4+N_DAMS] + w[5] - u[5+N_DAMS],
                x[6] - u[6] + u[2] + u[5] + u[2+N_DAMS] + u[5+N_DAMS] + w[6] - u[6+N_DAMS],
                x[7] - u[7] + w[7] - u[7+N_DAMS],
                x[8] - u[8] + u[6] + u[6+N_DAMS] + w[8] - u[8+N_DAMS]+ u[7] + u[7+N_DAMS],
                x[9] - u[9] + u[8] + u[8+N_DAMS] + w[9] - u[9+N_DAMS],
                x[10] - u[10] + u[9] + u[9+N_DAMS] + w[10] - u[10+N_DAMS]]
    end

    if LINEAR_COST
        function cost_t(t, x, u, w)
            return - COST[t] * sum(u[1:N_DAMS])
        end
    else
        function cost_t(t, x, u, w)
            return (- COST[t] * sum(u[1:N_DAMS]) + EPSILON_U*sum(u[1:N_DAMS]'*u[1:N_DAMS]) +
                        EPSILON_Z*sum((u[1:N_DAMS]+u[N_DAMS+1:2*N_DAMS])'*(u[1:N_DAMS]+u[N_DAMS+1:2*N_DAMS])))
        end
    end

    function final_cost_dams(model, m)
        alpha = JuMP.getvariable(m, :alpha)
        w = JuMP.getvariable(m, :w)
        x = JuMP.getvariable(m, :x)
        u = JuMP.getvariable(m, :u)
        xf = JuMP.getvariable(m, :xf)
        @JuMP.variable(m, z1[1:N_DAMS] >= 0)
        @JuMP.constraint(m, alpha == 0.)
        @JuMP.constraint(m, z1 .>= STOCK_TARGET - xf)
        @JuMP.objective(m, Min, model.costFunctions(model.stageNumber-1, x, u, w) +
                                alpha + COST_COEFF*sum(z1'*z1))
    end

#------------------------------
# 12 dams valley
#------------------------------
elseif VALLEY == "damXII"
    const N_DAMS = 12
    # Damsvalley configuration:
    const STOCK_TARGET = [40 for i in 1:N_DAMS]
    const LOWER_DAM = [2, 6, 5, 5, 6, 8, 9, 11, 11, 12, 12, 13]

    const X_BOUNDS = [(0, 80),
                      (0, 100),
                      (0, 80),
                      (0, 80),
                      (0, 100),
                      (0, 100),
                      (0, 80),
                      (0, 100),
                      (0, 100),
                      (0, 80),
                      (0, 100),
                      (0, 100)];
    const U_BOUNDS = vcat([(0, 40),
                           (0, 40),
                           (0, 40),
                           (0, 40),
                           (0, 80),
                           (0, 120),
                           (0, 40),
                           (0, 120),
                           (0, 40),
                           (0, 40),
                           (0, 140),
                           (0, 140)],
                            [(0, Inf) for i in 1:N_DAMS]);

    # Discretization of controls for assessment:
    DU = [10, 10, 10, 10, 10, 20, 10, 20, 10, 10, 20, 20, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8]

    function dynamic(t, x, u, w)
        return [x[1] - u[1] + w[1] - u[1+N_DAMS],
                x[2] - u[2] + u[1] + u[1+N_DAMS] + w[2] - u[2+N_DAMS],
                x[3] - u[3] + w[3] - u[3+N_DAMS],
                x[4] - u[4] + w[4] - u[4+N_DAMS],
                x[5] - u[5] + u[3] + u[4] + u[3+N_DAMS] + u[4+N_DAMS] + w[5] - u[5+N_DAMS],
                x[6] - u[6] + u[2] + u[5] + u[2+N_DAMS] + u[5+N_DAMS] + w[6] - u[6+N_DAMS],
                x[7] - u[7] + w[7] - u[7+N_DAMS],
                x[8] - u[8] + u[6] + u[6+N_DAMS] + w[8] - u[8+N_DAMS],
                x[9] - u[9] + u[7] + u[7+N_DAMS] + w[9] - u[9+N_DAMS],
                x[10] - u[10] + w[10] - u[10+N_DAMS],
                x[11] - u[11] + u[8] + u[8+N_DAMS]+  u[9] + u[9+N_DAMS] + w[11] - u[11+N_DAMS],
                x[12] - u[12] + u[10] + u[10+N_DAMS]+  u[11] + u[11+N_DAMS] + w[12] - u[12+N_DAMS]]
    end

    if LINEAR_COST
        function cost_t(t, x, u, w)
            return - COST[t] * sum(u[1:N_DAMS])
        end
    else
        function cost_t(t, x, u, w)
            return (- COST[t] * sum(u[1:N_DAMS]) + EPSILON_U*sum(u[1:N_DAMS]'*u[1:N_DAMS]) +
                        EPSILON_Z*sum((u[1:N_DAMS]+u[N_DAMS+1:2*N_DAMS])'*(u[1:N_DAMS]+u[N_DAMS+1:2*N_DAMS])))
        end
    end

    function final_cost_dams(model, m)
        alpha = JuMP.getvariable(m, :alpha)
        w = JuMP.getvariable(m, :w)
        x = JuMP.getvariable(m, :x)
        u = JuMP.getvariable(m, :u)
        xf = JuMP.getvariable(m, :xf)
        @JuMP.variable(m, z1[1:N_DAMS] >= 0)
        @JuMP.constraint(m, alpha == 0.)
        @JuMP.constraint(m, z1 .>= STOCK_TARGET - xf)
        @JuMP.objective(m, Min, model.costFunctions(model.stageNumber-1, x, u, w) +
                                alpha + COST_COEFF*sum(z1'*z1))
    end

#------------------------------
# 14 dams valley
#------------------------------
elseif VALLEY == "damXIV"
    const N_DAMS = 14
    # Damsvalley configuration:
    const STOCK_TARGET = [40 for i in 1:N_DAMS]
    const LOWER_DAM = [2, 6, 5, 5, 6, 8, 9, 11, 11, 12, 12, 13, 14, 15]

    const X_BOUNDS = [(0, 80),
                      (0, 100),
                      (0, 80),
                      (0, 80),
                      (0, 100),
                      (0, 100),
                      (0, 80),
                      (0, 100),
                      (0, 100),
                      (0, 80),
                      (0, 100),
                      (0, 100),
                      (0, 100),
                      (0, 100)];
    const U_BOUNDS = vcat([(0, 40),
                           (0, 40),
                           (0, 40),
                           (0, 40),
                           (0, 80),
                           (0, 120),
                           (0, 40),
                           (0, 120),
                           (0, 40),
                           (0, 40),
                           (0, 200),
                           (0, 220),
                           (0, 240),
                           (0, 260)],
                            [(0, Inf) for i in 1:N_DAMS]);

    # Discretization of controls for assessment:
    DU = [10, 10, 10, 10, 10, 20, 10, 20, 10, 10, 20, 20, 20, 20, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8]

    function dynamic(t, x, u, w)
        return [x[1] - u[1] + w[1] - u[1+N_DAMS],
                x[2] - u[2] + u[1] + u[1+N_DAMS] + w[2] - u[2+N_DAMS],
                x[3] - u[3] + w[3] - u[3+N_DAMS],
                x[4] - u[4] + w[4] - u[4+N_DAMS],
                x[5] - u[5] + u[3] + u[4] + u[3+N_DAMS] + u[4+N_DAMS] + w[5] - u[5+N_DAMS],
                x[6] - u[6] + u[2] + u[5] + u[2+N_DAMS] + u[5+N_DAMS] + w[6] - u[6+N_DAMS],
                x[7] - u[7] + w[7] - u[7+N_DAMS],
                x[8] - u[8] + u[6] + u[6+N_DAMS] + w[8] - u[8+N_DAMS],
                x[9] - u[9] + u[7] + u[7+N_DAMS] + w[9] - u[9+N_DAMS],
                x[10] - u[10] + w[10] - u[10+N_DAMS],
                x[11] - u[11] + u[8] + u[8+N_DAMS]+  u[9] + u[9+N_DAMS] + w[11] - u[11+N_DAMS],
                x[12] - u[12] + u[10] + u[10+N_DAMS]+  u[11] + u[11+N_DAMS] + w[12] - u[12+N_DAMS],
                x[13] - u[13] + u[12] + u[12+N_DAMS] + w[13] - u[13+N_DAMS],
                x[14] - u[14] + u[13] + u[13+N_DAMS] + w[14] - u[14+N_DAMS]]
    end

    if LINEAR_COST
        function cost_t(t, x, u, w)
            return - COST[t] * sum(u[1:N_DAMS])
        end
    else
        function cost_t(t, x, u, w)
            return (- COST[t] * sum(u[1:N_DAMS]) + EPSILON_U*sum(u[1:N_DAMS]'*u[1:N_DAMS]) +
                        EPSILON_Z*sum((u[1:N_DAMS]+u[N_DAMS+1:2*N_DAMS])'*(u[1:N_DAMS]+u[N_DAMS+1:2*N_DAMS])))
        end
    end

    function final_cost_dams(model, m)
        alpha = JuMP.getvariable(m, :alpha)
        w = JuMP.getvariable(m, :w)
        x = JuMP.getvariable(m, :x)
        u = JuMP.getvariable(m, :u)
        xf = JuMP.getvariable(m, :xf)
        @JuMP.variable(m, z1[1:N_DAMS] >= 0)
        @JuMP.constraint(m, alpha == 0.)
        @JuMP.constraint(m, z1 .>= STOCK_TARGET - xf)
        @JuMP.objective(m, Min, model.costFunctions(model.stageNumber-1, x, u, w) +
                                alpha + COST_COEFF*sum(z1'*z1))
    end

#------------------------------
# 16 dams valley
#------------------------------
elseif VALLEY == "damXVI"
    const N_DAMS = 16
    # Damsvalley configuration:
    const STOCK_TARGET = [40 for i in 1:N_DAMS]
    const LOWER_DAM = [2, 6, 5, 5, 6, 8, 9, 11, 11, 12, 12, 14, 14, 15, 16, 17]

    const X_BOUNDS = [(0, 80),
                      (0, 100),
                      (0, 80),
                      (0, 80),
                      (0, 100),
                      (0, 100),
                      (0, 80),
                      (0, 100),
                      (0, 100),
                      (0, 80),
                      (0, 100),
                      (0, 100),
                      (0, 80),
                      (0, 100),
                      (0, 100),
                      (0, 100)];
    const U_BOUNDS = vcat([(0, 40),
                           (0, 40),
                           (0, 40),
                           (0, 40),
                           (0, 80),
                           (0, 120),
                           (0, 40),
                           (0, 120),
                           (0, 40),
                           (0, 40),
                           (0, 200),
                           (0, 220),
                           (0, 40),
                           (0, 240),
                           (0, 260),
                           (0, 260)],
                            [(0, Inf) for i in 1:N_DAMS]);

    # Discretization of controls for assessment:
    DU = [10, 10, 10, 10, 10, 20, 10, 20, 10, 10, 20, 20, 10, 20, 20, 20, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8]

    function dynamic(t, x, u, w)
        return [x[1] - u[1] + w[1] - u[1+N_DAMS],
                x[2] - u[2] + u[1] + u[1+N_DAMS] + w[2] - u[2+N_DAMS],
                x[3] - u[3] + w[3] - u[3+N_DAMS],
                x[4] - u[4] + w[4] - u[4+N_DAMS],
                x[5] - u[5] + u[3] + u[4] + u[3+N_DAMS] + u[4+N_DAMS] + w[5] - u[5+N_DAMS],
                x[6] - u[6] + u[2] + u[5] + u[2+N_DAMS] + u[5+N_DAMS] + w[6] - u[6+N_DAMS],
                x[7] - u[7] + w[7] - u[7+N_DAMS],
                x[8] - u[8] + u[6] + u[6+N_DAMS] + w[8] - u[8+N_DAMS],
                x[9] - u[9] + u[7] + u[7+N_DAMS] + w[9] - u[9+N_DAMS],
                x[10] - u[10] + w[10] - u[10+N_DAMS],
                x[11] - u[11] + u[8] + u[8+N_DAMS]+  u[9] + u[9+N_DAMS] + w[11] - u[11+N_DAMS],
                x[12] - u[12] + u[10] + u[10+N_DAMS]+  u[11] + u[11+N_DAMS] + w[12] - u[12+N_DAMS],
                x[13] - u[13] + w[13] - u[13+N_DAMS],
                x[14] - u[14] + u[12] + u[12+N_DAMS] + u[13] + u[13+N_DAMS] + w[14] - u[14+N_DAMS],
                x[15] - u[15] + u[14] + u[14+N_DAMS] + w[15] - u[15+N_DAMS],
                x[16] - u[16] + u[15] + u[15+N_DAMS] + w[16] - u[16+N_DAMS]]
    end

    if LINEAR_COST
        function cost_t(t, x, u, w)
            return - COST[t] * sum(u[1:N_DAMS])
        end
    else
        function cost_t(t, x, u, w)
            return (- COST[t] * sum(u[1:N_DAMS]) + EPSILON_U*sum(u[1:N_DAMS]'*u[1:N_DAMS]) +
                        EPSILON_Z*sum((u[1:N_DAMS]+u[N_DAMS+1:2*N_DAMS])'*(u[1:N_DAMS]+u[N_DAMS+1:2*N_DAMS])))
        end
    end

    function final_cost_dams(model, m)
        alpha = JuMP.getvariable(m, :alpha)
        w = JuMP.getvariable(m, :w)
        x = JuMP.getvariable(m, :x)
        u = JuMP.getvariable(m, :u)
        xf = JuMP.getvariable(m, :xf)
        @JuMP.variable(m, z1[1:N_DAMS] >= 0)
        @JuMP.constraint(m, alpha == 0.)
        @JuMP.constraint(m, z1 .>= STOCK_TARGET - xf)
        @JuMP.objective(m, Min, model.costFunctions(model.stageNumber-1, x, u, w) +
                                alpha + COST_COEFF*sum(z1'*z1))
    end

#------------------------------
# 18 dams valley
#------------------------------
elseif VALLEY == "damXVIII"
    const N_DAMS = 18
    # Damsvalley configuration:
    const STOCK_TARGET = [40 for i in 1:N_DAMS]
    const LOWER_DAM = [2, 6, 5, 5, 6, 8, 9, 11, 11, 12, 12, 14, 14, 16, 17, 17, 18, 19]

    const X_BOUNDS = [(0, 80),
                      (0, 100),
                      (0, 80),
                      (0, 80),
                      (0, 100),
                      (0, 100),
                      (0, 80),
                      (0, 100),
                      (0, 100),
                      (0, 80),
                      (0, 100),
                      (0, 100),
                      (0, 80),
                      (0, 100),
                      (0, 80),
                      (0, 100),
                      (0, 100),
                      (0, 100)];
    const U_BOUNDS = vcat([(0, 40),
                           (0, 40),
                           (0, 40),
                           (0, 40),
                           (0, 80),
                           (0, 120),
                           (0, 40),
                           (0, 120),
                           (0, 40),
                           (0, 40),
                           (0, 200),
                           (0, 220),
                           (0, 40),
                           (0, 240),
                           (0, 40),
                           (0, 260),
                           (0, 280),
                           (0, 280)],
                            [(0, Inf) for i in 1:N_DAMS]);

    # Discretization of controls for assessment:
    DU = [10, 10, 10, 10, 10, 20, 10, 20, 10, 10, 20, 20, 10, 20, 10, 20, 20, 20, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8]

    function dynamic(t, x, u, w)
        return [x[1] - u[1] + w[1] - u[1+N_DAMS],
                x[2] - u[2] + u[1] + u[1+N_DAMS] + w[2] - u[2+N_DAMS],
                x[3] - u[3] + w[3] - u[3+N_DAMS],
                x[4] - u[4] + w[4] - u[4+N_DAMS],
                x[5] - u[5] + u[3] + u[4] + u[3+N_DAMS] + u[4+N_DAMS] + w[5] - u[5+N_DAMS],
                x[6] - u[6] + u[2] + u[5] + u[2+N_DAMS] + u[5+N_DAMS] + w[6] - u[6+N_DAMS],
                x[7] - u[7] + w[7] - u[7+N_DAMS],
                x[8] - u[8] + u[6] + u[6+N_DAMS] + w[8] - u[8+N_DAMS],
                x[9] - u[9] + u[7] + u[7+N_DAMS] + w[9] - u[9+N_DAMS],
                x[10] - u[10] + w[10] - u[10+N_DAMS],
                x[11] - u[11] + u[8] + u[8+N_DAMS]+  u[9] + u[9+N_DAMS] + w[11] - u[11+N_DAMS],
                x[12] - u[12] + u[10] + u[10+N_DAMS]+  u[11] + u[11+N_DAMS] + w[12] - u[12+N_DAMS],
                x[13] - u[13] + w[13] - u[13+N_DAMS],
                x[14] - u[14] + u[12] + u[12+N_DAMS] + u[13] + u[13+N_DAMS] + w[14] - u[14+N_DAMS],
                x[15] - u[15] + w[15] - u[15+N_DAMS],
                x[16] - u[16] + u[14] + u[14+N_DAMS] + w[16] - u[16+N_DAMS],
                x[17] - u[17] + u[15] + u[15+N_DAMS] + u[16] + u[16+N_DAMS] + w[17] - u[17+N_DAMS],
                x[18] - u[18] + u[17] + u[17+N_DAMS] + w[18] - u[18+N_DAMS]]
    end

    if LINEAR_COST
        function cost_t(t, x, u, w)
            return - COST[t] * sum(u[1:N_DAMS])
        end
    else
        function cost_t(t, x, u, w)
            return (- COST[t] * sum(u[1:N_DAMS]) + EPSILON_U*sum(u[1:N_DAMS]'*u[1:N_DAMS]) +
                        EPSILON_Z*sum((u[1:N_DAMS]+u[N_DAMS+1:2*N_DAMS])'*(u[1:N_DAMS]+u[N_DAMS+1:2*N_DAMS])))
        end
    end

    function final_cost_dams(model, m)
        alpha = JuMP.getvariable(m, :alpha)
        w = JuMP.getvariable(m, :w)
        x = JuMP.getvariable(m, :x)
        u = JuMP.getvariable(m, :u)
        xf = JuMP.getvariable(m, :xf)
        @JuMP.variable(m, z1[1:N_DAMS] >= 0)
        @JuMP.constraint(m, alpha == 0.)
        @JuMP.constraint(m, z1 .>= STOCK_TARGET - xf)
        @JuMP.objective(m, Min, model.costFunctions(model.stageNumber-1, x, u, w) +
                                alpha + COST_COEFF*sum(z1'*z1))
    end

#------------------------------
# 20 dams valley
#------------------------------
elseif VALLEY == "damXX"
    const N_DAMS = 20
    # Damsvalley configuration:
    const STOCK_TARGET = [40 for i in 1:N_DAMS]
    const LOWER_DAM = [2, 6, 5, 5, 6, 8, 9, 11, 11, 12, 12, 14, 14, 16, 17, 17, 18, 19, 20, 21]

    const X_BOUNDS = [(0, 80),
                      (0, 100),
                      (0, 80),
                      (0, 80),
                      (0, 100),
                      (0, 100),
                      (0, 80),
                      (0, 100),
                      (0, 100),
                      (0, 80),
                      (0, 100),
                      (0, 100),
                      (0, 80),
                      (0, 100),
                      (0, 80),
                      (0, 100),
                      (0, 100),
                      (0, 100),
                      (0, 100),
                      (0, 100)];
    const U_BOUNDS = vcat([(0, 40),
                           (0, 40),
                           (0, 40),
                           (0, 40),
                           (0, 80),
                           (0, 120),
                           (0, 40),
                           (0, 120),
                           (0, 40),
                           (0, 40),
                           (0, 200),
                           (0, 220),
                           (0, 40),
                           (0, 240),
                           (0, 40),
                           (0, 260),
                           (0, 280),
                           (0, 280),
                           (0, 300),
                           (0, 300)],
                            [(0, Inf) for i in 1:N_DAMS]);

    # Discretization of controls for assessment:
    DU = [10, 10, 10, 10, 10, 20, 10, 20, 10, 10, 20, 20, 10, 20, 10, 20, 20, 20, 20, 20, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8]

    function dynamic(t, x, u, w)
        return [x[1] - u[1] + w[1] - u[1+N_DAMS],
                x[2] - u[2] + u[1] + u[1+N_DAMS] + w[2] - u[2+N_DAMS],
                x[3] - u[3] + w[3] - u[3+N_DAMS],
                x[4] - u[4] + w[4] - u[4+N_DAMS],
                x[5] - u[5] + u[3] + u[4] + u[3+N_DAMS] + u[4+N_DAMS] + w[5] - u[5+N_DAMS],
                x[6] - u[6] + u[2] + u[5] + u[2+N_DAMS] + u[5+N_DAMS] + w[6] - u[6+N_DAMS],
                x[7] - u[7] + w[7] - u[7+N_DAMS],
                x[8] - u[8] + u[6] + u[6+N_DAMS] + w[8] - u[8+N_DAMS],
                x[9] - u[9] + u[7] + u[7+N_DAMS] + w[9] - u[9+N_DAMS],
                x[10] - u[10] + w[10] - u[10+N_DAMS],
                x[11] - u[11] + u[8] + u[8+N_DAMS]+  u[9] + u[9+N_DAMS] + w[11] - u[11+N_DAMS],
                x[12] - u[12] + u[10] + u[10+N_DAMS]+  u[11] + u[11+N_DAMS] + w[12] - u[12+N_DAMS],
                x[13] - u[13] + w[13] - u[13+N_DAMS],
                x[14] - u[14] + u[12] + u[12+N_DAMS] + u[13] + u[13+N_DAMS] + w[14] - u[14+N_DAMS],
                x[15] - u[15] + w[15] - u[15+N_DAMS],
                x[16] - u[16] + u[14] + u[14+N_DAMS] + w[16] - u[16+N_DAMS],
                x[17] - u[17] + u[15] + u[15+N_DAMS] + u[16] + u[16+N_DAMS] + w[17] - u[17+N_DAMS],
                x[18] - u[18] + u[17] + u[17+N_DAMS] + w[18] - u[18+N_DAMS],
                x[19] - u[19] + u[18] + u[18+N_DAMS] + w[19] - u[19+N_DAMS],
                x[20] - u[20] + u[19] + u[19+N_DAMS] + w[20] - u[20+N_DAMS]]
    end

    if LINEAR_COST
        function cost_t(t, x, u, w)
            return - COST[t] * sum(u[1:N_DAMS])
        end
    else
        function cost_t(t, x, u, w)
            return (- COST[t] * sum(u[1:N_DAMS]) + EPSILON_U*sum(u[1:N_DAMS]'*u[1:N_DAMS]) +
                        EPSILON_Z*sum((u[1:N_DAMS]+u[N_DAMS+1:2*N_DAMS])'*(u[1:N_DAMS]+u[N_DAMS+1:2*N_DAMS])))
        end
    end

    function final_cost_dams(model, m)
        alpha = JuMP.getvariable(m, :alpha)
        w = JuMP.getvariable(m, :w)
        x = JuMP.getvariable(m, :x)
        u = JuMP.getvariable(m, :u)
        xf = JuMP.getvariable(m, :xf)
        @JuMP.variable(m, z1[1:N_DAMS] >= 0)
        @JuMP.constraint(m, alpha == 0.)
        @JuMP.constraint(m, z1 .>= STOCK_TARGET - xf)
        @JuMP.objective(m, Min, model.costFunctions(model.stageNumber-1, x, u, w) +
                                alpha + COST_COEFF*sum(z1'*z1))
    end

#------------------------------
# 25 dams valley
#------------------------------
elseif VALLEY == "damXXV"
    const N_DAMS = 25
    # Damsvalley configuration:
    const STOCK_TARGET = [40 for i in 1:N_DAMS]
    const LOWER_DAM = [2, 6, 5, 5, 6, 8, 9, 11, 11, 12, 12, 14, 14, 16, 17, 17, 18, 19, 24, 21, 22, 23, 24, 25, 26]

    const X_BOUNDS = [(0, 80),
                      (0, 100),
                      (0, 80),
                      (0, 80),
                      (0, 100),
                      (0, 100),
                      (0, 80),
                      (0, 100),
                      (0, 100),
                      (0, 80),
                      (0, 100),
                      (0, 100),
                      (0, 80),
                      (0, 100),
                      (0, 80),
                      (0, 100),
                      (0, 100),
                      (0, 100),
                      (0, 100),
                      (0, 80),
                      (0, 100),
                      (0, 100),
                      (0, 100),
                      (0, 100),
                      (0, 100)];
    const U_BOUNDS = vcat([(0, 40),
                           (0, 40),
                           (0, 40),
                           (0, 40),
                           (0, 80),
                           (0, 120),
                           (0, 40),
                           (0, 120),
                           (0, 40),
                           (0, 40),
                           (0, 200),
                           (0, 220),
                           (0, 40),
                           (0, 240),
                           (0, 40),
                           (0, 260),
                           (0, 280),
                           (0, 280),
                           (0, 300),
                           (0, 40),
                           (0, 40),
                           (0, 80),
                           (0, 120),
                           (0, 400),
                           (0, 400)],
                            [(0, Inf) for i in 1:N_DAMS]);

    # Discretization of controls for assessment:
    DU = [10, 10, 10, 10, 10, 20, 10, 20, 10, 10, 20, 20, 10, 20, 10, 20, 20, 20, 20, 10, 10, 20, 20, 20, 20, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8]

    function dynamic(t, x, u, w)
        return [x[1] - u[1] + w[1] - u[1+N_DAMS],
                x[2] - u[2] + u[1] + u[1+N_DAMS] + w[2] - u[2+N_DAMS],
                x[3] - u[3] + w[3] - u[3+N_DAMS],
                x[4] - u[4] + w[4] - u[4+N_DAMS],
                x[5] - u[5] + u[3] + u[4] + u[3+N_DAMS] + u[4+N_DAMS] + w[5] - u[5+N_DAMS],
                x[6] - u[6] + u[2] + u[5] + u[2+N_DAMS] + u[5+N_DAMS] + w[6] - u[6+N_DAMS],
                x[7] - u[7] + w[7] - u[7+N_DAMS],
                x[8] - u[8] + u[6] + u[6+N_DAMS] + w[8] - u[8+N_DAMS],
                x[9] - u[9] + u[7] + u[7+N_DAMS] + w[9] - u[9+N_DAMS],
                x[10] - u[10] + w[10] - u[10+N_DAMS],
                x[11] - u[11] + u[8] + u[8+N_DAMS]+  u[9] + u[9+N_DAMS] + w[11] - u[11+N_DAMS],
                x[12] - u[12] + u[10] + u[10+N_DAMS]+  u[11] + u[11+N_DAMS] + w[12] - u[12+N_DAMS],
                x[13] - u[13] + w[13] - u[13+N_DAMS],
                x[14] - u[14] + u[12] + u[12+N_DAMS] + u[13] + u[13+N_DAMS] + w[14] - u[14+N_DAMS],
                x[15] - u[15] + w[15] - u[15+N_DAMS],
                x[16] - u[16] + u[14] + u[14+N_DAMS] + w[16] - u[16+N_DAMS],
                x[17] - u[17] + u[15] + u[15+N_DAMS] + u[16] + u[16+N_DAMS] + w[17] - u[17+N_DAMS],
                x[18] - u[18] + u[17] + u[17+N_DAMS] + w[18] - u[18+N_DAMS],
                x[19] - u[19] + u[18] + u[18+N_DAMS] + w[19] - u[19+N_DAMS],
                x[20] - u[20] + w[20] - u[20+N_DAMS],
                x[21] - u[21] + u[20] + u[20+N_DAMS] + w[21] - u[21+N_DAMS],
                x[22] - u[22] + u[21] + u[21+N_DAMS] + w[22] - u[22+N_DAMS],
                x[23] - u[23] + u[22] + u[22+N_DAMS] + w[23] - u[23+N_DAMS],
                x[24] - u[24] + u[19] + u[19+N_DAMS] + u[23] + u[23+N_DAMS] + w[24] - u[24+N_DAMS],
                x[25] - u[25] + u[24] + u[24+N_DAMS] + w[25] - u[25+N_DAMS]]
    end

    if LINEAR_COST
        function cost_t(t, x, u, w)
            return - COST[t] * sum(u[1:N_DAMS])
        end
    else
        function cost_t(t, x, u, w)
            return (- COST[t] * sum(u[1:N_DAMS]) + EPSILON_U*sum(u[1:N_DAMS]'*u[1:N_DAMS]) +
                        EPSILON_Z*sum((u[1:N_DAMS]+u[N_DAMS+1:2*N_DAMS])'*(u[1:N_DAMS]+u[N_DAMS+1:2*N_DAMS])))
        end
    end

    function final_cost_dams(model, m)
        alpha = JuMP.getvariable(m, :alpha)
        w = JuMP.getvariable(m, :w)
        x = JuMP.getvariable(m, :x)
        u = JuMP.getvariable(m, :u)
        xf = JuMP.getvariable(m, :xf)
        @JuMP.variable(m, z1[1:N_DAMS] >= 0)
        @JuMP.constraint(m, alpha == 0.)
        @JuMP.constraint(m, z1 .>= STOCK_TARGET - xf)
        @JuMP.objective(m, Min, model.costFunctions(model.stageNumber-1, x, u, w) +
                                alpha + COST_COEFF*sum(z1'*z1))
    end

#------------------------------
# Vicdessos Valley (DamRC1)
#------------------------------
elseif VALLEY == "damRC1"
    const N_DAMS = 5
    # Damsvalley configuration:
    const STOCK_TARGET = [15, 5, 15, 3, 3]
    const LOWER_DAM = [4, 4, 4, 5, 6]
    const X_BOUNDS = [(0, 30), (0, 10), (0, 30), (0, 6), (0, 6)];
    const U_BOUNDS = vcat([(0, 10), (0, 8), (0, 20), (0, 40), (0, 40)],
                            [(0., Inf) for i in 1:N_DAMS]);

    # Discretization of controls for assessment:
    DU = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2]

    function dynamic(t, x, u, w)
        return [x[1] - u[1] + w[1] - u[6],
                x[2] - u[2] + w[2] - u[7],
                x[3] - u[3] + w[3] - u[8],
                x[4] - u[4] + u[3] + u[2] + u[1] + u[8] + u[7] + u[6] + w[4] - u[9],
                x[5] - u[5] + u[4] + u[9] + w[5] - u[10]]
    end

    if LINEAR_COST
        function cost_t(t, x, u, w)
            return - COST[t] * (u[1] + u[2] + u[3] + u[4] + u[5])
        end
    else
        function cost_t(t, x, u, w)
            return - COST[t] * (u[1] + u[2] + u[3] + u[4] + u[5]) + EPSILON_U*(u[1]*u[1] + u[2]*u[2] + u[3]*u[3] + u[4]*u[4] + u[5]*u[5]) +
                        EPSILON_Z*((u[1]+u[6])*(u[1]+u[6]) + (u[2]+u[7])*(u[2]+u[7]) + (u[3]+u[8])*(u[3]+u[8])+ (u[4]+u[9])*(u[4]+u[9])+ (u[5]+u[10])*(u[5]+u[10]))
        end
    end

    function final_cost_dams(model, m)
        alpha = JuMP.getvariable(m, :alpha)
        w = JuMP.getvariable(m, :w)
        x = JuMP.getvariable(m, :x)
        u = JuMP.getvariable(m, :u)
        xf = JuMP.getvariable(m, :xf)
        @JuMP.variable(m, z1 >= 0)
        @JuMP.variable(m, z2 >= 0)
        @JuMP.variable(m, z3 >= 0)
        @JuMP.variable(m, z4 >= 0)
        @JuMP.variable(m, z5 >= 0)
        @JuMP.constraint(m, alpha == 0.)
        @JuMP.constraint(m, z1 >= STOCK_TARGET[1] - xf[1])
        @JuMP.constraint(m, z2 >= STOCK_TARGET[2] - xf[2])
        @JuMP.constraint(m, z3 >= STOCK_TARGET[3] - xf[3])
        @JuMP.constraint(m, z4 >= STOCK_TARGET[4] - xf[4])
        @JuMP.constraint(m, z5 >= STOCK_TARGET[5] - xf[5])
        @JuMP.objective(m, Min, model.costFunctions(model.stageNumber-1, x, u, w) + alpha + COST_COEFF*(z1*z1+z2*z2+z3*z3+z4*z4+z5*z5))
    end

#------------------------------
# Dordogne Valley (DamRC2)
#------------------------------
elseif VALLEY == "damRC2"
    const N_DAMS = 5
    # Damsvalley configuration:
    const STOCK_TARGET = [100, 20, 80, 60, 20]
    const LOWER_DAM = [2, 3, 4, 5, 6]
    const X_BOUNDS = [(0, 200), (0, 40), (0, 160), (0, 120), (0, 40)];
    const U_BOUNDS = vcat([(0, 400), (0, 500), (0, 600), (0, 700), (0, 800)],
                            [(0, Inf) for i in 1:N_DAMS]);

    # Discretization of controls for assessment:
    DU = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10]

    function dynamic(t, x, u, w)
        return [x[1] - u[1] + w[1] - u[6],
                x[2] - u[2] + u[1] + u[6] + w[2] - u[7],
                x[3] - u[3] + u[2] + u[7] + w[3] - u[8],
                x[4] - u[4] + u[3] + u[8] + w[4] - u[9],
                x[5] - u[5] + u[4] + u[9] + w[5] - u[10]]
    end

    if LINEAR_COST
        function cost_t(t, x, u, w)
            return - COST[t] * (u[1] + u[2] + u[3] + u[4] + u[5])
        end
    else
        function cost_t(t, x, u, w)
            return - COST[t] * (u[1] + u[2] + u[3] + u[4] + u[5]) + EPSILON_U*(u[1]*u[1] + u[2]*u[2] + u[3]*u[3] + u[4]*u[4] + u[5]*u[5]) +
                        EPSILON_Z*((u[1]+u[6])*(u[1]+u[6]) + (u[2]+u[7])*(u[2]+u[7]) + (u[3]+u[8])*(u[3]+u[8])+ (u[4]+u[9])*(u[4]+u[9])+ (u[5]+u[10])*(u[5]+u[10]))
        end
    end

    function final_cost_dams(model, m)
        alpha = JuMP.getvariable(m, :alpha)
        w = JuMP.getvariable(m, :w)
        x = JuMP.getvariable(m, :x)
        u = JuMP.getvariable(m, :u)
        xf = JuMP.getvariable(m, :xf)
        @JuMP.variable(m, z1 >= 0)
        @JuMP.variable(m, z2 >= 0)
        @JuMP.variable(m, z3 >= 0)
        @JuMP.variable(m, z4 >= 0)
        @JuMP.variable(m, z5 >= 0)
        @JuMP.constraint(m, alpha == 0.)
        @JuMP.constraint(m, z1 >= STOCK_TARGET[1] - xf[1])
        @JuMP.constraint(m, z2 >= STOCK_TARGET[2] - xf[2])
        @JuMP.constraint(m, z3 >= STOCK_TARGET[3] - xf[3])
        @JuMP.constraint(m, z4 >= STOCK_TARGET[4] - xf[4])
        @JuMP.constraint(m, z5 >= STOCK_TARGET[5] - xf[5])
        @JuMP.objective(m, Min, model.costFunctions(model.stageNumber-1, x, u, w) + alpha + COST_COEFF*(z1*z1+z2*z2+z3*z3+z4*z4+z5*z5))
    end

else
    println("Valley $VALLEY is not implemented")
end
