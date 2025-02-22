import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.StringTokenizer;

// Y = MT*D + max(B)*D
// MА = MT*(MT+MZ) - MZ*MT

public class MESIF_Variant1 {
    private double[][] MT;
    private double[][] MZ;
    private double[] B;
    private double[] D;

    public static void main(String[] args) {
        System.out.println("Var 1");

        MESIF_Variant1 var1 = new MESIF_Variant1();
        Map<Integer, Long> formula1Times = new LinkedHashMap<>();
        Map<Integer, Long> formula2Times = new LinkedHashMap<>();

        for (int size = 100; size <= 1000; size += 20) {
            final int currentSize = size;
            String filePath = "generated_data/data_" + size + ".txt";
            System.out.println("Reading file: " + filePath);
            var1.loadData(filePath, size);

            String resultFilePath = "results_var1/size_" + currentSize + ".txt";

            try (BufferedWriter writer = new BufferedWriter(new FileWriter(resultFilePath))) {
                writer.write("");
            } catch (IOException e) {
                System.err.println("Error clearing file: " + e.getMessage());
            }

            Thread formula1Thread = new Thread(() -> {
                double[] localD = var1.D.clone();
                double[] localB = var1.B.clone();
                double[][] localMT = var1.MT.clone();

                long startTime = System.nanoTime();
                double[] Y = calc_formula1(localMT, localD, localB);
                long endTime = System.nanoTime();
                formula1Times.put(currentSize, endTime - startTime);

                writeResultsToFile(Y, resultFilePath, "Y");
                printVector(Y);
            });

            Thread formula2Thread = new Thread(() -> {
                double[][] localMT = var1.MT.clone();
                double[][] localMZ = var1.MZ.clone();

                long startTime = System.nanoTime();
                double[][] MA = calc_formula2(localMT, localMZ);
                long endTime = System.nanoTime();
                formula2Times.put(currentSize, endTime - startTime);

                writeResultsToFile(MA, resultFilePath, "MA");
                printMatrix(MA);
            });

            formula1Thread.start();
            formula2Thread.start();

            try {
                formula1Thread.join();
                formula2Thread.join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        writeExecutionTimesToFile("results_var1/times.txt", formula1Times, formula2Times);

        System.out.println("Size,Formula1(ns),Formula2(ns)");
        for (Integer size : formula1Times.keySet()) {
            System.out.println(size + "," + formula1Times.get(size) + "," + formula2Times.get(size));
        }
    }

    public static void writeResultsToFile(double[] vector, String filePath, String name) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath, true))) {
            writer.write("Vector " + name + ":\n");
            for (double value : vector) {
                writer.write(value + "\n");
            }
            writer.write("\n");
        } catch (IOException e) {
            System.err.println("Error writing file: " + e.getMessage());
        }
    }

    public static void writeResultsToFile(double[][] matrix, String filePath, String name) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath, true))) {
            writer.write("Matrix " + name + ":\n");
            for (double[] row : matrix) {
                for (double value : row) {
                    writer.write(value + " ");
                }
                writer.newLine();
            }
            writer.write("\n");
        } catch (IOException e) {
            System.err.println("Error writing file: " + e.getMessage());
        }
    }

    public static void writeExecutionTimesToFile(String filePath, Map<Integer, Long> formula1Times,
            Map<Integer, Long> formula2Times) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            writer.write("Size,Formula1(ns),Formula2(ns)\n");
            for (Integer size : formula1Times.keySet()) {
                writer.write(size + "," + formula1Times.get(size) + "," + formula2Times.get(size) + "\n");
            }
        } catch (IOException e) {
            System.err.println("Error writing execution times file: " + e.getMessage());
        }
    }

    public static void printVector(double[] vector) {
        System.out.println("Vector:");
        for (double value : vector) {
            System.out.print(value + " ");
        }
        System.out.println();
    }

    public static void printMatrix(double[][] matrix) {
        System.out.println("Matrix:");
        for (double[] row : matrix) {
            for (double value : row) {
                System.out.print(value + " ");
            }
            System.out.println();
        }
    }

    public void loadData(String filePath, int size) {
        MT = new double[size][size];
        MZ = new double[size][size];
        B = new double[size];
        D = new double[size];

        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = reader.readLine()) != null) {
                StringTokenizer tokenizer = new StringTokenizer(line);
                if (!tokenizer.hasMoreTokens())
                    continue;

                String key = tokenizer.nextToken();
                int rows = Integer.parseInt(tokenizer.nextToken());
                int cols = Integer.parseInt(tokenizer.nextToken());

                if (key.equals("MT")) {
                    readMatrix(reader, MT, rows, cols);
                } else if (key.equals("MZ")) {
                    readMatrix(reader, MZ, rows, cols);
                } else if (key.equals("B")) {
                    readVector(reader, B, rows);
                } else if (key.equals("D")) {
                    readVector(reader, D, rows);
                }
            }
        } catch (IOException e) {
            System.err.println("Error reading file: " + e.getMessage());
        }
    }

    private static void readMatrix(BufferedReader reader, double[][] matrix, int rows, int cols) throws IOException {
        for (int i = 0; i < rows; i++) {
            StringTokenizer tokenizer = new StringTokenizer(reader.readLine());
            for (int j = 0; j < cols; j++) {
                matrix[i][j] = Double.parseDouble(tokenizer.nextToken().replace(",", "."));
            }
        }
    }

    private static void readVector(BufferedReader reader, double[] vector, int size) throws IOException {
        for (int i = 0; i < size; i++) {
            vector[i] = Double.parseDouble(reader.readLine().trim().replace(",", "."));
        }
    }

    private static double[] multiplyMatrixVectorKahan(double[][] matrix, double[] vector) {
        int size = vector.length;
        double[] result = new double[size];
        for (int i = 0; i < size; i++) {
            double sum = 0.0;
            double c = 0.0;
            for (int j = 0; j < size; j++) {
                double y = matrix[i][j] * vector[j] - c;
                double t = sum + y;
                c = (t - sum) - y;
                sum = t;
            }
            result[i] = sum;
        }
        return result;
    }

    private static double[][] multiplyMatricesKahan(double[][] matrix1, double[][] matrix2) {
        int size = matrix1.length;
        double[][] result = new double[size][size];
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                double sum = 0.0;
                double c = 0.0;
                for (int k = 0; k < size; k++) {
                    double y = matrix1[i][k] * matrix2[k][j] - c;
                    double t = sum + y;
                    c = (t - sum) - y;
                    sum = t;
                }
                result[i][j] = sum;
            }
        }
        return result;
    }

    private static double kahanSum(double a, double b) {
        double sum = 0.0;
        double c = 0.0;
        double y = b - c;
        double t = a + y;
        c = (t - a) - y;
        sum = t;
        return sum;
    }

    private static double[][] addMatricesKahan(double[][] matrix1, double[][] matrix2) {
        int size = matrix1.length;
        double[][] result = new double[size][size];
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                result[i][j] = kahanSum(matrix1[i][j], matrix2[i][j]);
            }
        }
        return result;
    }

    private static double findMax(double[] vector) {
        return Arrays.stream(vector).max().orElse(0);
    }

    public static double[] calc_formula1(double[][] MT_local, double[] D_local, double[] B_local) {
        double[] MT_D = multiplyMatrixVectorKahan(MT_local, D_local);
        double maxB = findMax(B_local);
        double[] Y = new double[D_local.length];
        for (int i = 0; i < D_local.length; i++) {
            Y[i] = kahanSum(MT_D[i], maxB * D_local[i]);
        }
        return Y;
    }

    public static double[][] calc_formula2(double[][] MT_local, double[][] MZ_local) {
        double[][] MT_plus_MZ = addMatricesKahan(MT_local, MZ_local);
        double[][] MT_MT_plus_MZ = multiplyMatricesKahan(MT_local, MT_plus_MZ);
        double[][] MZ_MT = multiplyMatricesKahan(MZ_local, MT_local);
        int size = MT_local.length;
        double[][] MA = new double[size][size];
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                MA[i][j] = kahanSum(MT_MT_plus_MZ[i][j], -MZ_MT[i][j]);
            }
        }
        return MA;
    }
}
