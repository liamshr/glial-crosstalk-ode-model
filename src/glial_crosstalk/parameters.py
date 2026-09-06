"""Model parameters and initial conditions taken from the original notebook."""

from dataclasses import asdict, dataclass, replace


@dataclass(frozen=True)
class ModelParameters:
    k_n1: float = 1.23121468628907
    k_n12: float = 1.3444424885940314
    x_n: float = 2.47585316397624
    b_n: float = 0.8682932268535773
    a_inf1: float = 2.944948653936038
    k_n4: float = 1.272493506699988
    k_n10: float = 0.8577380244036067
    z_n: float = 2.3786552539888537
    y_n: float = 0.9996433532137823
    s_mr: float = 1.0960359353301747
    mu_mr: float = 0.36232256985896183
    tau_n4: float = 1.1705406202322937
    tau_n10: float = 0.9048457485245999
    g_n: float = 2.464228895378905
    m_n: float = 0.563495322717254
    k_M1base: float = 0.9926354332956264
    kcd: float = 1.202623143443761
    h_n: float = 2.493063875169106
    v_n: float = 1.6633392059433056
    a_inf2: float = 1.651861798448556
    k_pn1: float = 2.00487443066468
    mu_n1: float = 0.31766812384668947
    k_pn12: float = 20.659169915759403
    mu_n12: float = 0.20508819119658298
    k_tbase: float = 0.947187881552978
    k_tn12: float = 1.0707050213093836
    c_n: float = 2.661117852785608
    r_n: float = 1.6774534754063852
    k_M2base: float = 0.9707774602596897
    kc4: float = 0.9238073930606325
    q_n: float = 2.3567649771794734
    w_n: float = 0.9997073991548225
    k_tn10: float = 1.2260936900650916
    k_pn10: float = 1.2436568657494924
    mu_n10: float = 0.47020795161687023
    k_tn4: float = 0.9546207188594067
    k_pn4: float = 0.9803768557089533
    mu_n4: float = 0.30439335856121924
    alpha_n12: float = 1.1383684645435397
    alpha_n1: float = 1.1915264432034482
    r_M1: float = 0.26862120846959575
    gamma_M1: float = 0.38999592743353073
    gamma_M2: float = 0.5308181970140535
    mu_M1: float = 0.09685420836766151
    mu_M2: float = 0.08041242746608893
    k_Aq_M2: float = 0.3
    k_Ap_M1: float = 0.8
    k_Ap_D: float = 0.6
    k_Aq2Ap: float = 0.4
    k_Ap2Aq: float = 0.2
    mu_Aq: float = 0.05
    mu_Ap: float = 0.08
    s_Aq: float = 0.1
    beta_Aq_M2: float = 0.5
    beta_Ap_M1: float = 0.4
    beta_Aq_inhib_M1: float = 0.3
    psi_Aq: float = 0.3
    psi_Ap: float = 0.2
    ROS: float = 100.0
    p_ros: float = 1.0
    i_ros: float = 10.0

    def updated(self, **changes):
        """Return a validated immutable parameter set with selected changes."""
        return replace(self, **changes)

    def as_dict(self):
        return asdict(self)


DEFAULT_PARAMETERS = ModelParameters()

# (M1, M2, IL1, IL12, IL10, IL4, Aq, Ap, D)
DEFAULT_INITIAL_CONDITIONS = (0.1, 0.1, 0.1, 0.01, 0.1, 0.01, 0.05, 0.005, 0.0)
