package au.org.ala.spatial.analysis.layers;

import org.apache.commons.math3.util.Pair;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Map;

public class ShannonHLayerGenerator extends CalculatedLayerGenerator {

    public ShannonHLayerGenerator(BigDecimal resolution, File coordinateSpeciesFlatFile) throws IOException {
        super(resolution);
        readCoordinateSpeciesFlatFile(coordinateSpeciesFlatFile);
    }

    public static void main(String[] args) throws Exception {
        if (args.length == 0) {
            System.out.println("args[0]=Resolution in degrees from the list; 1, 0.1, 0.01\n"
                    + "args[1]=Path to coordinate,species flat file (should be generated by biocache store via jenkins - file resolution must match the resolution provided to this tool)\n"
                    + "args[2]=Path of directory in which to write output files\n"
                    + "args[3]=Prefix to use for names of output files.\n");

            return;
        }

        BigDecimal resolution = new BigDecimal(args[0]).setScale(2);
        File coordinateSpeciesFlatFile = new File(args[1]);
        File outputFileDirectory = new File(args[2]);
        String outputFileNamePrefix = args[3];

        new ShannonHLayerGenerator(resolution, coordinateSpeciesFlatFile).writeGrid(outputFileDirectory, outputFileNamePrefix);
    }

    @Override
    protected float handleCell(Pair<BigDecimal, BigDecimal> coordPair, float maxValue, PrintWriter ascPrintWriter, BufferedOutputStream divaOutputStream) throws IOException {
        if (_cellSpeciesOccurrenceCounts.containsKey(coordPair)) {
            // Calculate shannon H value for the cell. -1 * Sum[i=1-n] (pi * ln (pi))
            // where i is each species in the cell
            // and pi is number of occurrences of species i / total number of occurrences in the cell

            Map<Integer, Integer> m = _cellSpeciesOccurrenceCounts.get(coordPair);

            int totalOccurrences = 0;
            for (Integer i : m.values()) {
                totalOccurrences += i;
            }

            double sum = 0;

            for (Integer i : m.values()) {
                double p = i / (double) totalOccurrences;
                sum += p * Math.log(p);
            }

            sum *= -1;

            ascPrintWriter.print((float) sum);

            ByteBuffer bb = ByteBuffer.wrap(new byte[Float.SIZE / Byte.SIZE]);
            bb.order(ByteOrder.LITTLE_ENDIAN);
            bb.putFloat((float) sum);
            divaOutputStream.write(bb.array());

            return maxValue < sum ? (float) sum : maxValue;
        } else {
            // No species occurrences in this cell.
            // is zero.
            ascPrintWriter.print("0");

            ByteBuffer bb = ByteBuffer.wrap(new byte[Float.SIZE / Byte.SIZE]);
            bb.order(ByteOrder.LITTLE_ENDIAN);
            bb.putFloat(0);
            divaOutputStream.write(bb.array());

            return maxValue;
        }
    }
}
