#include "generalizedcubic.hpp"
#include <cmath>

//  **********  This code is preliminary, and will be updated.  **********
//  **********  Use only for beta testing.  **********

double AbstractCubic::am_term(double tau, std::size_t itau)
{
    double summer = 0;
    for(int i = N-1; i >= 0; --i)
    {
        for (int j = N-1; j >= 0; --j)
        {
            summer += x[i]*x[j]*aij_term(tau, i, j, itau);
        }
    }
    return summer;
}
double AbstractCubic::d_am_term_dxi(double tau, std::size_t itau, std::size_t i, bool xN_independent)
{
    if (xN_independent)
    {
        double summer = 0;
        for (int j = N-1; j >= 0; --j)
        {
            summer += x[j]*aij_term(tau, i, j, itau);
        }
        return 2*summer;
    }
    else{
        double summer = 0;
        for (int k = N-2; k >= 0; --k)
        {
            summer += x[k]*(aij_term(tau, i, k, itau)-aij_term(tau, k, N-1, itau));
        }
        return 2*(summer + x[N-1]*(aij_term(tau, N-1, i, itau) - aij_term(tau, N-1, N-1, itau)));
    }
}
double AbstractCubic::d2_am_term_dxidxj(double tau, std::size_t itau, std::size_t i, std::size_t j, bool xN_independent)
{
    if (xN_independent)
    {
        return 2*aij_term(tau, i, j, itau);
    }
    else{
        return 2*(aij_term(tau, i, j, itau)-aij_term(tau, j, N-1, itau)-aij_term(tau, N-1, i, itau)+aij_term(tau, N-1, N-1, itau));
    }
}

double AbstractCubic::d3_am_term_dxidxjdxk(double tau, std::size_t itau, std::size_t i, std::size_t j, std::size_t k, bool xN_independent)
{
    return 0;
}

double AbstractCubic::bm_term()
{
    double summer = 0;
    for(int i = N-1; i >= 0; --i)
    {
        summer += x[i]*b0_ii(i);
    }
    return summer;
}
double AbstractCubic::d_bm_term_dxi(std::size_t i, bool xN_independent)
{
    if (xN_independent)
    {
        return b0_ii(i);
    }
    else{
        return b0_ii(i) - b0_ii(N-1);
    }
}
double AbstractCubic::d2_bm_term_dxidxj(std::size_t i, std::size_t j, bool xN_independent)
{
    return 0;
}
double AbstractCubic::d3_bm_term_dxidxjdxk(std::size_t i, std::size_t j, std::size_t k, bool xN_independent)
{
    return 0;
}

double AbstractCubic::aii_term(double tau, std::size_t i, std::size_t itau)
{
    double Tr_over_Tci = T_r/Tc[i];
    double sqrt_Tr_Tci = sqrt(Tr_over_Tci);
    // If we are not using the full Mathias-Copeman formulation for a_ii,
    // we just use the simple results from the supplemental information because
    // they are much more computationally efficient
    if (simple_aii){
        // All derivatives have a common bracketed term, so we factor it out
        // and calculate it here
        double m = m_ii(i);
        double B = 1 + m*(1-sqrt_Tr_Tci*sqrt(1/tau));

        switch (itau){
        case 0:
            return a0_ii(i)*B*B;
        case 1:
            return a0_ii(i)*m*B/pow(tau, 3.0/2.0)*sqrt_Tr_Tci;
        case 2:
            return a0_ii(i)*m/2.0*(m/pow(tau, 3)*Tr_over_Tci - 3*B/pow(tau, 5.0/2.0)*sqrt_Tr_Tci);
        case 3:
            return (3.0/4.0)*a0_ii(i)*m*(-3.0*m/pow(tau, 4)*Tr_over_Tci + 5*B/pow(tau, 7.0/2.0)*sqrt_Tr_Tci);
        case 4:
            return (3.0/8.0)*a0_ii(i)*m*(29.0*m/pow(tau, 5)*Tr_over_Tci - 35*B/pow(tau, 9.0/2.0)*sqrt_Tr_Tci);
        default:
            throw -1;
        }
    }
    else{
        // Here we are using the full Mathias-Copeman formulation, introducing
        // some additional computational effort, so we only evaluate the parameters that
        // we actually need to evaluate, otherwise we just set their value to zero
        // See info on the conditional (ternary) operator : http://www.cplusplus.com/articles/1AUq5Di1/
        // Furthermore, this should help with branch prediction
        double Di = 1-sqrt_Tr_Tci/sqrt(tau);
        double dDi_dtau = (itau >= 1) ? (1.0/2.0)*sqrt_Tr_Tci/(pow(tau,1.5)) : 0;
        double d2Di_dtau2 = (itau >= 2) ? -(3.0/4.0)*sqrt_Tr_Tci/(pow(tau,2.5)) : 0;
        double d3Di_dtau3 = (itau >= 3) ? (15.0/8.0)*sqrt_Tr_Tci/(pow(tau,3.5)) : 0;
        double d4Di_dtau4 = (itau >= 4) ? -(105.0/16.0)*sqrt_Tr_Tci/(pow(tau,4.5)) : 0;

        double Bi = 1, dBi_dtau = 0, d2Bi_dtau2 = 0, d3Bi_dtau3 = 0, d4Bi_dtau4 = 0;
        for (int n = 1; n <= 3; ++n){
            const std::vector<double> &C = get_C_ref(static_cast<int>(n));
            Bi += C[i]*pow(Di, n);
            dBi_dtau += (itau < 1) ? 0 : (n*C[i]*pow(Di,n-1)*dDi_dtau) ;
            d2Bi_dtau2 += (itau < 2) ? 0 : n*C[i]*((n-1)*pow(dDi_dtau,2) + Di*d2Di_dtau2)*pow(Di, n-2);
            d3Bi_dtau3 += (itau < 3) ? 0 : n*C[i]*(3*(n-1)*Di*dDi_dtau*d2Di_dtau2 + (n*n-3*n+2)*pow(dDi_dtau,3)+pow(Di,2)*d3Di_dtau3)*pow(Di, n-3);
            d4Bi_dtau4 += (itau < 4) ? 0 : n*C[i]*(6*(n*n-3*n+2)*Di*pow(dDi_dtau,2)*d2Di_dtau2 + (n*n*n-6*n*n+11*n-6)*pow(dDi_dtau,4)
                                                   +(4*n*dDi_dtau*d3Di_dtau3+3*n*pow(d2Di_dtau2,2)-4*dDi_dtau*d3Di_dtau3-3*pow(d2Di_dtau2,2) )*pow(Di,2)
                                                   + pow(Di,3)*d4Di_dtau4 )*pow(Di, n-4);
        }
        switch (itau){
        case 0:
            return a0_ii(i)*Bi*Bi;
        case 1:
            return 2*a0_ii(i)*Bi*dBi_dtau;
        case 2:
            return 2*a0_ii(i)*(Bi*d2Bi_dtau2 + dBi_dtau*dBi_dtau);
        case 3:
            return 2*a0_ii(i)*(Bi*d3Bi_dtau3 + 3*dBi_dtau*d2Bi_dtau2);
        case 4:
            return 2*a0_ii(i)*(Bi*d4Bi_dtau4 + 4*dBi_dtau*d3Bi_dtau3 + 3*pow(d2Bi_dtau2, 2));
        default:
            throw -1;
        }
    }
}
double AbstractCubic::u_term(double tau, std::size_t i, std::size_t j, std::size_t itau)
{
    double aii = aii_term(tau, i, 0), ajj = aii_term(tau, j, 0);
    switch (itau){
    case 0:
        return aii*ajj;
    case 1:
        return aii*aii_term(tau, j, 1) + ajj*aii_term(tau, i, 1);
    case 2:
        return (aii*aii_term(tau, j, 2)
               +2*aii_term(tau, i, 1)*aii_term(tau, j, 1)
               +ajj*aii_term(tau, i, 2)
               );
    case 3:
        return (aii*aii_term(tau, j, 3)
                +3*aii_term(tau, i, 1)*aii_term(tau, j, 2)
                +3*aii_term(tau, i, 2)*aii_term(tau, j, 1)
                +ajj*aii_term(tau, i, 3)
                );
    case 4:
        return (aii*aii_term(tau, j, 4)
                +4*aii_term(tau, i, 1)*aii_term(tau, j, 3)
                +6*aii_term(tau, i, 2)*aii_term(tau, j, 2)
                +4*aii_term(tau, i, 3)*aii_term(tau, j, 1)
                +ajj*aii_term(tau, i, 4)
                );
    default:
        throw -1;
    }
}
double AbstractCubic::aij_term(double tau, std::size_t i, std::size_t j, std::size_t itau)
{
    double u = u_term(tau, i, j, 0);

    switch (itau){
    case 0:
        return (1-k[i][j])*sqrt(u);
    case 1:
        return (1-k[i][j])/(2.0*sqrt(u))*u_term(tau, i, j, 1);
    case 2:
        return (1-k[i][j])/(4.0*pow(u,3.0/2.0))*(2*u*u_term(tau, i, j, 2)-pow(u_term(tau, i, j, 1), 2));
    case 3:
        return (1-k[i][j])/(8.0*pow(u,5.0/2.0))*(4*pow(u,2)*u_term(tau, i, j, 3)
                                                -6*u*u_term(tau, i, j, 1)*u_term(tau, i, j, 2)
                                                +3*pow(u_term(tau, i, j, 1),3));
    case 4:
        return (1-k[i][j])/(16.0*pow(u,7.0/2.0))*(-4*pow(u,2)*(4*u_term(tau, i, j, 1)*u_term(tau, i, j, 3) + 3*pow(u_term(tau, i, j, 2),2))
                                                +8*pow(u,3)*u_term(tau, i, j, 4) + 36*u*pow(u_term(tau, i, j, 1),2)*u_term(tau, i, j, 2)
                                                -15*pow(u_term(tau, i, j, 1), 4)
                                                );
    default:
        throw -1;
    }
}
double AbstractCubic::psi_minus(double delta, std::size_t itau, std::size_t idelta)
{
    if (itau > 0) return 0.0;
    double b = bm_term();
    double bracket = 1-b*delta*rho_r;

    switch(idelta){
    case 0:
        return -log(bracket);
    case 1:
        return b*rho_r/bracket;
    case 2:
        return pow(b*rho_r/bracket, 2);
    case 3:
        return 2*pow(b*rho_r/bracket, 3);
    case 4:
        return 6*pow(b*rho_r/bracket, 4);
    default:
        throw -1;
    }
}
double AbstractCubic::d_psi_minus_dxi(double delta, std::size_t itau, std::size_t idelta, std::size_t i, bool xN_independent)
{
    if (itau > 0) return 0.0;
    double b = bm_term();
    double db_dxi = d_bm_term_dxi(i, xN_independent);
    double bracket = 1-b*delta*rho_r;

    switch(idelta){
    case 0:
        return delta*rho_r*db_dxi/bracket;
    case 1:
        return rho_r*db_dxi/pow(bracket, 2);
    case 2:
        return 2*pow(rho_r,2)*b*db_dxi/pow(bracket, 3);
    case 3:
        return 6*pow(rho_r,3)*pow(b, 2)*db_dxi/pow(bracket, 4);
    case 4:
        return 24*pow(rho_r,4)*pow(b, 3)*db_dxi/pow(bracket, 5);
    default:
        throw -1;
    }
}
double AbstractCubic::d2_psi_minus_dxidxj(double delta, std::size_t itau, std::size_t idelta, std::size_t i, std::size_t j, bool xN_independent)
{
    if (itau > 0) return 0.0;
    double b = bm_term();
    double db_dxi = d_bm_term_dxi(i, xN_independent),
           db_dxj = d_bm_term_dxi(j, xN_independent),
           d2b_dxidxj = d2_bm_term_dxidxj(i, j, xN_independent);
    double bracket = 1-b*delta*rho_r;

    switch(idelta){
    case 0:
        return pow(delta*rho_r, 2)*db_dxi*db_dxj/pow(bracket, 2) + delta*rho_r*d2b_dxidxj/bracket;
    case 1:
        return 2*delta*pow(rho_r, 2)*db_dxi*db_dxj/pow(bracket, 3) + rho_r*d2b_dxidxj/pow(bracket, 2);
    case 2:
        return 2*pow(rho_r,2)*db_dxi*db_dxj/pow(bracket, 4)*(2*delta*rho_r*b+1) + 2*pow(rho_r, 2)*b*d2b_dxidxj/pow(bracket,3);
    case 3:
        return 12*pow(rho_r,3)*b*db_dxi*db_dxj/pow(bracket, 5)*(delta*rho_r*b+1) + 6*pow(rho_r, 3)*pow(b,2)*d2b_dxidxj/pow(bracket,4);
    case 4:
        return 24*pow(rho_r,4)*pow(b, 2)*db_dxi*db_dxj/pow(bracket, 6)*(2*delta*rho_r*b + 3) + 24*pow(rho_r, 4)*pow(b,3)*d2b_dxidxj/pow(bracket,5);
    default:
        throw -1;
    }
}
double AbstractCubic::d3_psi_minus_dxidxjdxk(double delta, std::size_t itau, std::size_t idelta, std::size_t i, std::size_t j, std::size_t k, bool xN_independent)
{
    if (itau > 0) return 0.0;
    double b = bm_term();
    double db_dxi = d_bm_term_dxi(i, xN_independent),
           db_dxj = d_bm_term_dxi(j, xN_independent),
           db_dxk = d_bm_term_dxi(k, xN_independent),
           d2b_dxidxj = d2_bm_term_dxidxj(i, j, xN_independent),
           d2b_dxidxk = d2_bm_term_dxidxj(i, k, xN_independent),
           d2b_dxjdxk = d2_bm_term_dxidxj(j, k, xN_independent),
           d3b_dxidxjdxk = d3_bm_term_dxidxjdxk(i, j, k, xN_independent);
    double bracket = 1-b*delta*rho_r;

    switch(idelta){
    case 0:
        return delta*rho_r*d3b_dxidxjdxk/bracket
               + 2*pow(delta*rho_r, 3)*db_dxi*db_dxj*db_dxk/pow(bracket, 3)
               + pow(delta*rho_r, 2)/pow(bracket, 2)*(db_dxi*d2b_dxjdxk
                                                      +db_dxj*d2b_dxidxk
                                                      +db_dxk*d2b_dxidxj);
    case 1:
        return rho_r*d3b_dxidxjdxk/pow(bracket, 2)
               + 6*pow(delta, 2)*pow(rho_r, 3)*db_dxi*db_dxj*db_dxk/pow(bracket, 4)
               + 2*delta*pow(rho_r, 2)/pow(bracket, 3)*(db_dxi*d2b_dxjdxk
                                                        +db_dxj*d2b_dxidxk
                                                        +db_dxk*d2b_dxidxj);
    default:
        throw -1;
    }
}
double AbstractCubic::PI_12(double delta, std::size_t idelta)
{
    double b = bm_term();
    switch(idelta){
    case 0:
        return (1+Delta_1*b*rho_r*delta)*(1+Delta_2*b*rho_r*delta);
    case 1:
        return b*rho_r*(2*Delta_1*Delta_2*b*delta*rho_r+Delta_1+Delta_2);
    case 2:
        return 2*Delta_1*Delta_2*pow(b*rho_r, 2);
    case 3:
        return 0;
    case 4:
        return 0;
    default:
        throw -1;
    }
}
double AbstractCubic::d_PI_12_dxi(double delta, std::size_t idelta, std::size_t i, bool xN_independent)
{
    double b = bm_term();
    double db_dxi = d_bm_term_dxi(i, xN_independent);
    switch(idelta){
    case 0:
        return delta*rho_r*db_dxi*(2*Delta_1*Delta_2*b*delta*rho_r+Delta_1+Delta_2);
    case 1:
        return rho_r*db_dxi*(4*Delta_1*Delta_2*b*delta*rho_r+Delta_1+Delta_2);
    case 2:
        return 4*Delta_1*Delta_2*pow(rho_r, 2)*b*db_dxi;
    case 3:
        return 0;
    case 4:
        return 0;
    default:
        throw -1;
    }
}
double AbstractCubic::d2_PI_12_dxidxj(double delta, std::size_t idelta, std::size_t i, std::size_t j, bool xN_independent)
{
    double b = bm_term();
    double db_dxi = d_bm_term_dxi(i, xN_independent),
           db_dxj = d_bm_term_dxi(j, xN_independent),
           d2b_dxidxj = d2_bm_term_dxidxj(i, j, xN_independent);
    switch(idelta){
    case 0:
        return delta*rho_r*(2*Delta_1*Delta_2*delta*rho_r*db_dxi*db_dxj + (2*Delta_1*Delta_2*delta*rho_r*b+Delta_1+Delta_2)*d2b_dxidxj);
    case 1:
        return rho_r*(4*Delta_1*Delta_2*delta*rho_r*db_dxi*db_dxj + (4*Delta_1*Delta_2*delta*rho_r*b+Delta_1+Delta_2)*d2b_dxidxj);
    case 2:
        return 4*Delta_1*Delta_2*pow(rho_r,2)*(db_dxi*db_dxj + b*d2b_dxidxj);
    case 3:
        return 0;
    case 4:
        return 0;
    default:
        throw -1;
    }
}
double AbstractCubic::d3_PI_12_dxidxjdxk(double delta, std::size_t idelta, std::size_t i, std::size_t j, std::size_t k, bool xN_independent)
{
    double b = bm_term();
    double db_dxi = d_bm_term_dxi(i, xN_independent),
           db_dxj = d_bm_term_dxi(j, xN_independent),
           db_dxk = d_bm_term_dxi(k, xN_independent),
           d2b_dxidxj = d2_bm_term_dxidxj(i, j, xN_independent),
           d2b_dxidxk = d2_bm_term_dxidxj(i, k, xN_independent),
           d2b_dxjdxk = d2_bm_term_dxidxj(j, k, xN_independent),
           d3b_dxidxjdxk = d3_bm_term_dxidxjdxk(i, j, k, xN_independent);
    switch(idelta){
    case 0:
        return delta*rho_r*((2*Delta_1*Delta_2*delta*rho_r*b+Delta_1+Delta_2)*d3b_dxidxjdxk
                            + 2*Delta_1*Delta_2*delta*rho_r*(db_dxi*d2b_dxjdxk
                            +db_dxj*d2b_dxidxk
                            +db_dxk*d2b_dxidxj
                            )
                            );
    case 1:
        return rho_r*((4*Delta_1*Delta_2*delta*rho_r*b+Delta_1+Delta_2)*d3b_dxidxjdxk
                       + 4*Delta_1*Delta_2*delta*rho_r*(db_dxi*d2b_dxjdxk
                       + db_dxj*d2b_dxidxk
                       + db_dxk*d2b_dxidxj
                        )
                    );
    default:
        throw -1;
    }
}
double AbstractCubic::psi_plus(double delta, std::size_t idelta)
{
    switch(idelta){
    case 0:
        return A_term(delta)*c_term()/(Delta_1-Delta_2);
    case 1:
        return rho_r/PI_12(delta,0);
    case 2:
        return -rho_r/pow(PI_12(delta,0),2)*PI_12(delta,1);
    case 3:
        return rho_r*(-PI_12(delta,0)*PI_12(delta,2)+2*pow(PI_12(delta,1),2))/pow(PI_12(delta,0),3);
    case 4:
        // Term -PI_12(delta,x,0)*PI_12(delta,x,3) in the numerator is zero (and removed) since PI_12(delta,x,3) = 0
        return rho_r*(6*PI_12(delta,0)*PI_12(delta,1)*PI_12(delta,2) - 6*pow(PI_12(delta,1),3))/pow(PI_12(delta,0),4);
    default:
        throw -1;
    }
}
double AbstractCubic::d_psi_plus_dxi(double delta, std::size_t idelta, std::size_t i, bool xN_independent)
{
    double bracket = 0;
    if (idelta == 0){
        return (A_term(delta)*d_c_term_dxi(i, xN_independent) + c_term()*d_A_term_dxi(delta, i, xN_independent))/(Delta_1 - Delta_2);
    }
    // All the terms with at least one delta derivative are multiplied by a common term of -rhor/PI12^2
    // So we just evaluated the bracketed term and then multiply by the common factor in the front
    switch(idelta){
    case 1:
        bracket = d_PI_12_dxi(delta, 0, i, xN_independent); break;
    case 2:
        bracket = (d_PI_12_dxi(delta, 1, i, xN_independent)
                   + 2/rho_r*PI_12(delta,0)*PI_12(delta,1)*d_psi_plus_dxi(delta, 1,i,xN_independent)
                   );
        break;
    case 3:{
        bracket = (d_PI_12_dxi(delta, 2, i, xN_independent)
                   + 2/rho_r*(pow(PI_12(delta, 1), 2) + PI_12(delta, 0)*PI_12(delta, 2))*d_psi_plus_dxi(delta, 1,i,xN_independent)
                   + 4/rho_r*PI_12(delta,0)*PI_12(delta,1)*d_psi_plus_dxi(delta,2,i,xN_independent)
                   );
        break;
    }
    case 4:
        // d_PI_12_dxi(delta, x, 3, i, xN_independent) = 0, and PI_12(delta,x,0)*PI_12(delta,x,3) = 0, so removed from sum
        bracket = (6/rho_r*PI_12(delta,1)*PI_12(delta,2)*d_psi_plus_dxi(delta,1,i,xN_independent)
                   + 6/rho_r*(pow(PI_12(delta,1), 2) + PI_12(delta,0)*PI_12(delta,2))*d_psi_plus_dxi(delta,2,i,xN_independent)
                   + 6/rho_r*PI_12(delta,0)*PI_12(delta,1)*d_psi_plus_dxi(delta,3,i,xN_independent)
                   );
        break;
    default:
        throw -1;
    }
    return -rho_r/pow(PI_12(delta,0), 2)*bracket;
}
double AbstractCubic::d2_psi_plus_dxidxj(double delta, std::size_t idelta, std::size_t i, std::size_t j, bool xN_independent)
{
    double bracket = 0;
    double PI12 = PI_12(delta, 0);
    if (idelta == 0){
        return (A_term(delta)*d2_c_term_dxidxj(i, j, xN_independent)
                +c_term()*d2_A_term_dxidxj(delta, i, j, xN_independent)
                +d_A_term_dxi(delta, i, xN_independent)*d_c_term_dxi(j, xN_independent)
                +d_A_term_dxi(delta, j, xN_independent)*d_c_term_dxi(i, xN_independent)
                )/(Delta_1 - Delta_2);
    }
    // All the terms with at least one delta derivative have a common factor of -1/PI_12^2 out front
    // so we just calculate the bracketed term and then multiply later on
    switch(idelta){
    case 1:
        bracket = (rho_r*d2_PI_12_dxidxj(delta, 0, i, j, xN_independent)
                   + 2*PI12*d_PI_12_dxi(delta, 0, j, xN_independent)*d_psi_plus_dxi(delta, 1, i, xN_independent)
                   );
        break;
    case 2:
        bracket = (rho_r*d2_PI_12_dxidxj(delta, 1, i, j, xN_independent)
                   + 2*(PI12*d_PI_12_dxi(delta, 1, j, xN_independent)
                        + PI_12(delta, 1)*d_PI_12_dxi(delta, 0, j, xN_independent)
                        )*d_psi_plus_dxi(delta, 1, i, xN_independent)
                   + 2*PI12*PI_12(delta, 1)*d2_psi_plus_dxidxj(delta, 1, i, j, xN_independent)
                   + 2*PI12*d_PI_12_dxi(delta, 0, j, xN_independent)*d_psi_plus_dxi(delta, 2, i, xN_independent)
                   );
        break;
    case 3:{
        bracket = (rho_r*d2_PI_12_dxidxj(delta, 2, i, j, xN_independent)
                   + 2*(PI12*PI_12(delta, 2) + pow(PI_12(delta, 1), 2))*d2_psi_plus_dxidxj(delta, 1, i, j, xN_independent)
                   + 4*(PI12*d_PI_12_dxi(delta, 1, j, xN_independent)
                        + PI_12(delta, 1)*d_PI_12_dxi(delta, 0, j, xN_independent)
                        )*d_psi_plus_dxi(delta, 2, i, xN_independent)
                   + 2*(  PI12*d_PI_12_dxi(delta, 2, j, xN_independent)
                          + 2*PI_12(delta, 1)*d_PI_12_dxi(delta, 1, j, xN_independent)
                          + d_PI_12_dxi(delta, 0, j, xN_independent)*PI_12(delta, 2)
                          )*d_psi_plus_dxi(delta, 1, i, xN_independent)
                   + 4*PI12*PI_12(delta, 1)*d2_psi_plus_dxidxj(delta, 2, i, j, xN_independent)
                   + 2*PI12*d_PI_12_dxi(delta, 0, j, xN_independent)*d_psi_plus_dxi(delta, 3, i, xN_independent)
                   );
        break;
    }
    case 4:
        // rho_r*d2_PI_12_dxidxj(delta, x, 3, i, j, xN_independent)  = 0
        // PI_12(delta, x, 3) = 0
        // PI12*d_PI_12_dxi(delta, x, 3, j, xN_independent) = 0
        // d_PI_12_dxi(delta, x, 0, j, xN_independent)*PI_12(delta, x, 3) = 0
        bracket = (
                   + 6*(PI12*PI_12(delta, 2) + pow(PI_12(delta, 1), 2))*d2_psi_plus_dxidxj(delta, 2, i, j, xN_independent)
                   + 6*PI_12(delta, 1)*PI_12(delta, 2)*d2_psi_plus_dxidxj(delta, 1, i, j, xN_independent)
                   + 6*(PI12*d_PI_12_dxi(delta, 1, j, xN_independent)
                        + PI_12(delta, 1)*d_PI_12_dxi(delta, 0, j, xN_independent)
                        )*d_psi_plus_dxi(delta, 3, i, xN_independent)
                   + 6*(PI12*d_PI_12_dxi(delta, 2, j, xN_independent)
                        + 2*PI_12(delta, 1)*d_PI_12_dxi(delta, 1, j, xN_independent)
                        + d_PI_12_dxi(delta, 0, j, xN_independent)*PI_12(delta, 2)
                        )*d_psi_plus_dxi(delta, 2, i, xN_independent)
                   + 6*(PI_12(delta, 1)*d_PI_12_dxi(delta, 2, j, xN_independent)
                        + PI_12(delta, 2)*d_PI_12_dxi(delta, 1, j, xN_independent)
                        )*d_psi_plus_dxi(delta, 1, i, xN_independent)
                   + 6*PI12*PI_12(delta, 1)*d2_psi_plus_dxidxj(delta, 3, i, j, xN_independent)
                   + 2*PI12*d_PI_12_dxi(delta, 0, j, xN_independent)*d_psi_plus_dxi(delta, 4, i, xN_independent)
                   );
        break;
    default:
        throw -1;
    }
    return -1/pow(PI12, 2)*bracket;
}
double AbstractCubic::d3_psi_plus_dxidxjdxk(double delta, std::size_t idelta, std::size_t i, std::size_t j, std::size_t k, bool xN_independent)
{
    double PI12 = PI_12(delta, 0);
    switch (idelta){
    case 0:
        return (A_term(delta)*d3_c_term_dxidxjdxk(i, j, k, xN_independent)
            +c_term()*d3_A_term_dxidxjdxk(delta, i, j, k, xN_independent)
            +d_A_term_dxi(delta, i, xN_independent)*d2_c_term_dxidxj(j, k, xN_independent)
            +d_A_term_dxi(delta, j, xN_independent)*d2_c_term_dxidxj(i, k, xN_independent)
            +d_A_term_dxi(delta, k, xN_independent)*d2_c_term_dxidxj(i, j, xN_independent)
            +d_c_term_dxi(i, xN_independent)*d2_A_term_dxidxj(delta, j, k, xN_independent)
            +d_c_term_dxi(j, xN_independent)*d2_A_term_dxidxj(delta, i, k, xN_independent)
            +d_c_term_dxi(k, xN_independent)*d2_A_term_dxidxj(delta, i, j, xN_independent)
            )/(Delta_1 - Delta_2);
    case 1:
        return -1/pow(PI12, 2)*(rho_r*d3_PI_12_dxidxjdxk(delta, 0, i, j, k, xN_independent)
                                +2*(PI12*d2_PI_12_dxidxj(delta, 0, j, k, xN_independent) + d_PI_12_dxi(delta, 0, j, xN_independent)*d_PI_12_dxi(delta, 0, k, xN_independent))*d_psi_plus_dxi(delta, 1, i, xN_independent)
                                +2*PI12*d_PI_12_dxi(delta, 0, j, xN_independent)*d2_psi_plus_dxidxj(delta, 1, i, k, xN_independent) + 2*PI12*d_PI_12_dxi(delta, 0, k, xN_independent)*d2_psi_plus_dxidxj(delta, 1, i, j, xN_independent)
                                );
    default:
        throw -1;
    }
}

double AbstractCubic::tau_times_a(double tau, std::size_t itau)
{
    if (itau == 0){
        return tau*am_term(tau, 0);
    }
    else{
        return tau*am_term(tau, itau) + itau*am_term(tau, itau-1);
    }
}
double AbstractCubic::d_tau_times_a_dxi(double tau, std::size_t itau, std::size_t i, bool xN_independent)
{
    if (itau == 0){
        return tau*d_am_term_dxi(tau, 0,i,xN_independent);
    }
    else{
        return tau*d_am_term_dxi(tau, itau,i,xN_independent) + itau*d_am_term_dxi(tau, itau-1,i,xN_independent);
    }
}
double AbstractCubic::d2_tau_times_a_dxidxj(double tau, std::size_t itau, std::size_t i, std::size_t j, bool xN_independent)
{
    if (itau == 0){
        return tau*d2_am_term_dxidxj(tau, 0,i,j,xN_independent);
    }
    else{
        return tau*d2_am_term_dxidxj(tau, itau,i,j,xN_independent) + itau*d2_am_term_dxidxj(tau, itau-1,i,j,xN_independent);
    }
}
double AbstractCubic::d3_tau_times_a_dxidxjdxk(double tau, std::size_t itau, std::size_t i, std::size_t j, std::size_t k, bool xN_independent)
{
    if (itau == 0){
        return tau*d3_am_term_dxidxjdxk(tau, 0,i,j,k,xN_independent);
    }
    else{
        return tau*d3_am_term_dxidxjdxk(tau, itau,i,j,k,xN_independent) + itau*d3_am_term_dxidxjdxk(tau, itau-1,i,j,k,xN_independent);
    }
}
double AbstractCubic::alphar(double tau, double delta, std::size_t itau, std::size_t idelta)
{
    return psi_minus(delta, itau, idelta)-tau_times_a(tau,itau)/(R_u*T_r)*psi_plus(delta,idelta);
}
double AbstractCubic::d_alphar_dxi(double tau, double delta, std::size_t itau, std::size_t idelta, std::size_t i, bool xN_independent)
{
    return (d_psi_minus_dxi(delta, itau, idelta, i, xN_independent)
            -1/(R_u*T_r)*(d_tau_times_a_dxi(tau, itau,i,xN_independent)*psi_plus(delta, idelta)
                          +tau_times_a(tau, itau)*d_psi_plus_dxi(delta, idelta,i,xN_independent)
                          )
            );
}
double AbstractCubic::d2_alphar_dxidxj(double tau, double delta, std::size_t itau, std::size_t idelta, std::size_t i, std::size_t j,bool xN_independent)
{
    return (d2_psi_minus_dxidxj(delta, itau, idelta, i, j, xN_independent)
            -1/(R_u*T_r)*(d2_tau_times_a_dxidxj(tau, itau,i,j,xN_independent)*psi_plus(delta, idelta)
                          +d_tau_times_a_dxi(tau, itau,i,xN_independent)*d_psi_plus_dxi(delta, idelta,j,xN_independent)
                          +d_tau_times_a_dxi(tau, itau,j,xN_independent)*d_psi_plus_dxi(delta, idelta,i,xN_independent)
                          +tau_times_a(tau, itau)*d2_psi_plus_dxidxj(delta, idelta,i,j,xN_independent)
                          )
            );
}
double AbstractCubic::d3_alphar_dxidxjdxk(double tau, double delta, std::size_t itau, std::size_t idelta, std::size_t i, std::size_t j, std::size_t k, bool xN_independent)
{
    return (d3_psi_minus_dxidxjdxk(delta, itau, idelta, i, j, k, xN_independent)
            -1/(R_u*T_r)*(d2_tau_times_a_dxidxj(tau,itau,i,j,xN_independent)*d_psi_plus_dxi(delta,idelta,k,xN_independent)
                          +d3_tau_times_a_dxidxjdxk(tau,itau,i,j,k,xN_independent)*psi_plus(delta,idelta)

                          +d_tau_times_a_dxi(tau,itau,i,xN_independent)*d2_psi_plus_dxidxj(delta,idelta,j,k,xN_independent)
                          +d2_tau_times_a_dxidxj(tau,itau,i,k,xN_independent)*d_psi_plus_dxi(delta,idelta,j,xN_independent)

                          +d_tau_times_a_dxi(tau,itau,j,xN_independent)*d2_psi_plus_dxidxj(delta,idelta,i,k,xN_independent)
                          +d2_tau_times_a_dxidxj(tau,itau,j,k,xN_independent)*d_psi_plus_dxi(delta,idelta,i,xN_independent)

                          +tau_times_a(tau,itau)*d3_psi_plus_dxidxjdxk(delta,idelta,i,j,k,xN_independent)
                          +d_tau_times_a_dxi(tau,itau,k, xN_independent)*d2_psi_plus_dxidxj(delta,idelta,i,j,xN_independent)
                          )
            );
}



double SRK::a0_ii(std::size_t i)
{
    // Values from Soave, 1972 (Equilibium constants from a ..)
    double a = 0.42747*R_u*R_u*Tc[i]*Tc[i]/Pc[i];
    return a;
}
double SRK::b0_ii(std::size_t i)
{
    // Values from Soave, 1972 (Equilibium constants from a ..)
    double b = 0.08664*R_u*Tc[i]/Pc[i];
    return b;
}
double SRK::m_ii(std::size_t i)
{
    // Values from Soave, 1972 (Equilibium constants from a ..)
    double omega = acentric[i];
    double m = 0.480 + 1.574*omega - 0.176*omega*omega;
    return m;
}


double PengRobinson::a0_ii(std::size_t i)
{
    double a = 0.45724*R_u*R_u*Tc[i]*Tc[i]/Pc[i];
    return a;
}
double PengRobinson::b0_ii(std::size_t i)
{
    double b = 0.07780*R_u*Tc[i]/Pc[i];
    return b;
}
double PengRobinson::m_ii(std::size_t i)
{
double m;

    double omega = acentric[i];
    if (omega <= 0.491)
        m = 0.37464 + 1.54226*omega - 0.26992*omega*omega;
    else
        m = 0.379642 + 1.48503*omega -  0.164423*omega*omega + + 0.016666*omega*omega*omega;
    return m;
}
