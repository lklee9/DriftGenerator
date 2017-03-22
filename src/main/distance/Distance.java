package main.distance;

import com.yahoo.labs.samoa.instances.WekaToSamoaInstanceConverter;

/**
 * Created by loongkuan on 26/03/16.
 **/

public abstract class Distance {
    WekaToSamoaInstanceConverter wekaConverter = new WekaToSamoaInstanceConverter();
    public abstract double findDistance(double[] p, double[] q, double[] weights);
    public double findDistance(double[] p, double[] q) {
        double[] weights = new double[p.length];
        for (int i = 0; i < p.length; i++) {
            weights[i] = 1.0;
        }
        return this.findDistance(p, q, weights);
    }
}
