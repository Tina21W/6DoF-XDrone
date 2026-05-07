function prop_omega_cmd = prop_speed_controller(v_b, sim)
%PROP_SPEED_CONTROLLER Command relative prop angular rate from body x speed.
%   Positive u_error increases the magnitude of the nominal prop speed while
%   preserving the spin direction configured in the vehicle file.

    prop_omega_cmd = sim.prop.Prop.omega;

    if ~isfield(sim.prop.Prop, 'control')
        return
    end

    prop_control = sim.prop.Prop.control;
    u_des = prop_control.u_des;
    if isfield(sim.options, 'control') && isfield(sim.options.control, 'u_des')
        u_des = sim.options.control.u_des;
    end

    u_error = u_des - v_b(1);

    prop_omega_min = prop_control.omega_min;
    prop_omega_max = prop_control.omega_max;
    u_error_scale = prop_control.u_error_scale;
    prop_omega_delta_max = prop_control.omega_delta_max;

    prop_spin_sign = sign(sim.prop.Prop.omega);
    if prop_spin_sign == 0
        prop_spin_sign = 1;
    end

    prop_omega_delta = prop_omega_delta_max * tanh(u_error / u_error_scale);
    prop_omega_mag = abs(sim.prop.Prop.omega) + prop_omega_delta;
    prop_omega_mag = min(max(prop_omega_mag, prop_omega_min), prop_omega_max);
    prop_omega_cmd = prop_spin_sign * prop_omega_mag;

end
