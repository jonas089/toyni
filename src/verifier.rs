use crate::prover::StarkProof;

pub struct StarkVerifier;

impl StarkVerifier {
    pub fn verify(&self, proof: &StarkProof) -> bool {
        // todo: FRI folding consistency check at x
        /*
        let lhs = d_next[i];
        let rhs = d_prev[i] + beta * d_prev[i + half];
        assert_eq(lhs, rhs);
        */

        // spot check at z
        // public parameters:
        //  g from original domain (trace size evaluation domain)
        //  vanishing poly z from original domain
        //  trace poly commitments over shifted domain (t0, t1, t2)

        // prover inputs:
        //  q(z), q(x)
        /*assert_eq!(
            fibonacci_constraint(
                trace_poly.evaluate(g * g * z),
                trace_poly.evaluate(g * z),
                trace_poly.evaluate(z),
            ) * boundary_constraint_1(z, g, trace_len)
                * boundary_constraint_2(z, g, trace_len),
            q_poly.evaluate(z) * z_poly.evaluate(z)
        );*/
        true
    }
}
