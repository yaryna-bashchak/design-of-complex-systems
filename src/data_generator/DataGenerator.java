package data_generator;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;
import java.util.Map;
import java.util.List;
import java.util.HashMap;
import java.util.Arrays;

public class DataGenerator {
    private static final int START_SIZE = 100;
    private static final int END_SIZE = 1000;
    private static final int STEP = 20;
    private static final String DIRECTORY_NAME = "generated_data";

    public static void main(String[] args) {
        List<String> dataToGenerate = Arrays.asList("MT", "D", "B", "MZ");

        File directory = new File(DIRECTORY_NAME);
        if (!directory.exists()) {
            directory.mkdir();
        }

        for (int size = START_SIZE; size <= END_SIZE; size += STEP) {
            String fileName = DIRECTORY_NAME + "/data_" + size + ".txt";
            try (FileWriter writer = new FileWriter(fileName)) {
                generateData(writer, dataToGenerate, size);
                System.out.println("Data has been saved to: " + fileName);
            } catch (IOException e) {
                System.err.println("Error writing to file " + fileName + ": " + e.getMessage());
            }
        }
    }

    private static void generateData(FileWriter writer, List<String> dataToGenerate, int size) throws IOException {
        Random random = new Random();

        Map<String, Integer[]> dimensions = new HashMap<>();
        dimensions.put("MT", new Integer[]{size, size});
        dimensions.put("MZ", new Integer[]{size, size});
        dimensions.put("D", new Integer[]{size, 1});
        dimensions.put("B", new Integer[]{size, 1});

        for (String key : dataToGenerate) {
            if (dimensions.containsKey(key)) {
                Integer[] dim = dimensions.get(key);
                writer.write(key + " " + dim[0] + " " + dim[1] + "\n");

                for (int i = 0; i < dim[0]; i++) {
                    for (int j = 0; j < dim[1]; j++) {
                        writer.write(String.format("%.10f ", generatePositiveFloatingPointNumber(random)));
                    }
                    writer.write("\n");
                }
            } else {
                System.out.println("Warning: " + key + " is not defined in the data structure.");
            }
        }
    }

    private static double generatePositiveFloatingPointNumber(Random random) {
        double[] scales = {1e-5, 1e-2, 1.0, 1e2, 1e5};
        double scale = scales[random.nextInt(scales.length)];
        return (random.nextDouble() * scale) + 1e-10;
    }
}
