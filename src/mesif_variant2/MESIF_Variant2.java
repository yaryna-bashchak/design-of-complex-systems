import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.concurrent.locks.ReentrantLock;

// Y = MT*D + max(B)*D
// MА = MT*(MT+MZ) - MZ*MT

public class MESIF_Variant2 {
    private double[][] MT;
    private double[][] MZ;
    private double[] B;
    private double[] D;
    private double[] Y;
    private double[][] MA;

    private static double b_max = Double.NEGATIVE_INFINITY;
    private static final ReentrantLock lock = new ReentrantLock();

    public static void main(String[] args) {
        System.out.println("Var 2");

        MESIF_Variant2 var2 = new MESIF_Variant2();
        Map<Integer, Long> formula1Times = new LinkedHashMap<>();
        Map<Integer, Long> formula2Times = new LinkedHashMap<>();

        for (int size = 100; size <= 1000; size += 20) {
            final int currentSize = size;
            String filePath = "generated_data/data_" + size + ".txt";
            System.out.println("Reading file: " + filePath);
            var2.loadData(filePath, size);

            String resultFilePath = "results_var2/size_" + currentSize + ".txt";

            try (BufferedWriter writer = new BufferedWriter(new FileWriter(resultFilePath))) {
                writer.write("");
            } catch (IOException e) {
                System.err.println("Error clearing file: " + e.getMessage());
            }

            // formula 1
            long startTimeFormula1 = System.nanoTime();
            var2.calc_formula1();
            long endTimeFormula1 = System.nanoTime();
            formula1Times.put(currentSize, endTimeFormula1 - startTimeFormula1);

            writeResultsToFile(var2.Y, resultFilePath, "Y");
            printVector(var2.Y);

            // formula 2
            long startTimeFormula2 = System.nanoTime();
            var2.calc_formula2();
            long endTimeFormula2 = System.nanoTime();
            formula2Times.put(currentSize, endTimeFormula2 - startTimeFormula2);

            writeResultsToFile(var2.MA, resultFilePath, "MA");
            printMatrix(var2.MA);
        }

        writeExecutionTimesToFile("results_var2/times.txt", formula1Times, formula2Times);

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

    private static double kahanSum(double a, double b) {
        double sum = 0.0;
        double c = 0.0;
        double y = b - c;
        double t = a + y;
        c = (t - a) - y;
        sum = t;
        return sum;
    }

    private void multiplyPartMatrixVectorKahan(double[][] matrix, double[] vector, double[] result, int start, int end) {
        int size = vector.length;

        for (int i = start; i < end; i++) {
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
    }

    private void multiplyMatricesPartKahan(double[][] A, double[][] B, double[][] result, int start, int end) {
        for (int i = 0; i < B.length; i++) {
            for (int j = start; j < end; j++) {
                double sum = 0.0;
                double c = 0.0;
                for (int k = 0; k < B.length; k++) {
                    double y = A[i][k] * B[k][j] - c;
                    double t = sum + y;
                    c = (t - sum) - y;
                    sum = t;
                }
                result[i][j] = sum;
            }
        }
    }

    private void subtractMatricesPartKahan(double[][] A, double[][] B, double[][] result, int start, int end) {
        for (int i = 0; i < A.length; i++) {
            for (int j = start; j < end; j++) {
                result[i][j] = kahanSum(A[i][j], -B[i][j]);
            }
        }
    }

    public void calc_formula1() {
        int halfSize = B.length / 2;
        Y = new double[B.length];

        Thread t1 = new Thread(() -> {
            double localMax = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < halfSize; i++) {
                localMax = Math.max(localMax, B[i]);
            }
            lock.lock();
            try {
                b_max = Math.max(b_max, localMax);
            } finally {
                lock.unlock();
            }
        });

        Thread t2 = new Thread(() -> {
            double localMax = Double.NEGATIVE_INFINITY;
            for (int i = halfSize; i < B.length; i++) {
                localMax = Math.max(localMax, B[i]);
            }
            lock.lock();
            try {
                b_max = Math.max(b_max, localMax);
            } finally {
                lock.unlock();
            }
        });

        t1.start();
        t2.start();
        try {
            t1.join();
            t2.join();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        double[] MT_D = new double[B.length];

        t1 = new Thread(() -> {
            multiplyPartMatrixVectorKahan(MT, D, MT_D, 0, halfSize);
            for (int i = 0; i < halfSize; i++) {
                Y[i] = kahanSum(MT_D[i], b_max * D[i]);
            }
        });

        t2 = new Thread(() -> {
            multiplyPartMatrixVectorKahan(MT, D, MT_D, halfSize, B.length);
            for (int i = halfSize; i < B.length; i++) {
                Y[i] = kahanSum(MT_D[i], b_max * D[i]);
            }
        });

        t1.start();
        t2.start();
        try {
            t1.join();
            t2.join();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    public void calc_formula2() {
        int halfSize = MT.length / 2;
        double[][] temp1 = new double[MT.length][MT.length]; // (MT + MZ)
        double[][] temp2 = new double[MT.length][MT.length]; // MT * (MT + MZ)
        double[][] temp3 = new double[MT.length][MT.length]; // MZ * MT
        MA = new double[MT.length][MT.length];

        // temp1 = (MT + MZ)
        Thread t1 = new Thread(() -> {
            for (int i = 0; i < MT.length; i++) {
                for (int j = 0; j < halfSize; j++) {
                    temp1[i][j] = MT[i][j] + MZ[i][j]; // temp1 = (MT + MZ)
                }
            }
            multiplyMatricesPartKahan(MT, temp1, temp2, 0, halfSize); // temp2 = MT * temp1
            multiplyMatricesPartKahan(MZ, MT, temp3, 0, halfSize); // temp3 = MZ * MT
            subtractMatricesPartKahan(temp2, temp3, MA, 0, halfSize); // MA = temp2 - temp3
        });

        Thread t2 = new Thread(() -> {
            for (int i = 0; i < MT.length; i++) {
                for (int j = halfSize; j < MT.length; j++) {
                    temp1[i][j] = MT[i][j] + MZ[i][j];
                }
            }
            multiplyMatricesPartKahan(MT, temp1, temp2, halfSize, MT.length);
            multiplyMatricesPartKahan(MZ, MT, temp3, halfSize, MT.length);
            subtractMatricesPartKahan(temp2, temp3, MA, halfSize, MT.length);
        });

        t1.start();
        t2.start();
        try {
            t1.join();
            t2.join();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }
}
