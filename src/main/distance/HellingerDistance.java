package main.distance;

/**
 * Created by loongkuan on 26/03/16.
 **/
public class HellingerDistance extends Distance{
    @Override
    public double findDistance(double[] p, double[] q, double[] weights) {
        assert p.length == q.length;
        double driftMag = 0.0f;
        for (int i = 0; i < p.length; i++) {
            driftMag += Math.pow((Math.sqrt(p[i]) - Math.sqrt(q[i])), 2) * weights[i];
        }
        driftMag = Math.sqrt(driftMag);
        driftMag *= 1 / Math.sqrt(2);
        return driftMag;
    }
}
