%% @file f.m
%  @brief Computes the system of differential equations for an RLC circuit.
%  @author Hubert Durnas
%  
%  @param t Time variable [s]
%  @param Y State vector [u(t); v(t)]
%  @param R Resistance [Ω]
%  @param L Inductance [H]
%  @param C Capacitance [F]
%  @param U Input voltage [V]
%  @return y Vector containing first and second derivatives

function y = f(t, Y, R, L, C, U)
    % Initialize output vector
    y = zeros(2, 1);

    % @brief First equation: Voltage derivative
    y(1) = Y(2);

    % @brief Second equation: Current derivative
    y(2) = -R / L * Y(2) - 1 / (L * C) * Y(1) + U / (L * C);
end
