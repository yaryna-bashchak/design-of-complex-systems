package lab2_tsx;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

// Y = MT*D + max(B)*D
// MА = MT*(MT+MZ) - MZ*MT

public class Lab2 {
    private static final int NUM_THREADS = 2;
    private static double[][] MT;
    private static double[][] MZ;
    private static double[] B;
    private static double[] D;
    private static double[] Y;
    private static double[][] MA;
    private static double[] MT_D;
    private static double[][] MT_plus_MZ;
    private static double[][] MT_mlt_MT_plus_MZ;
    private static double[][] MZ_mlt_MT;
    private static double b_max = Double.NEGATIVE_INFINITY;

    private static final Lock lockY = new ReentrantLock();
    private static final Lock lockMA = new ReentrantLock();
    private static final Lock lockB_max = new ReentrantLock();
    private static final Lock lockMT_D = new ReentrantLock();
    private static final Lock lockMT_plus_MZ = new ReentrantLock();
    private static final Lock lockMT_mlt_MT_plus_MZ = new ReentrantLock();
    private static final Lock lockMZ_mlt_MT = new ReentrantLock();

    private static CyclicBarrier barrier;

    public static void main(String[] args) {
        Map<Integer, Long> times = new LinkedHashMap<>();

        for (int size = 100; size <= 1000; size += 20) {
            final int currentSize = size;
            barrier = new CyclicBarrier(NUM_THREADS);
            String filePath = "generated_data/data_" + size + ".txt";
            System.out.println("Reading file: " + filePath);
            loadData(filePath, size);

            String resultFilePath = "results_lab2_wo_tsx/size_" + currentSize + ".txt";

            try (BufferedWriter writer = new BufferedWriter(new FileWriter(resultFilePath))) {
                writer.write("");
            } catch (IOException e) {
                System.err.println("Error clearing file: " + e.getMessage());
            }

            Thread[] threads = new Thread[NUM_THREADS];
            int chunkSize = size / NUM_THREADS;
            
            long startTime = System.nanoTime();

            for (int i = 0; i < NUM_THREADS; i++) {
                int start = i * chunkSize;
                int end = (i == NUM_THREADS - 1) ? size - 1 : (start + chunkSize - 1);
                threads[i] = new Thread(new RunnableThread(start, end));
                threads[i].start();
            }
            
            for (Thread thread : threads) {
                try {
                    thread.join();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
            
            long endTime = System.nanoTime();
            times.put(currentSize, endTime - startTime);

            writeResultsToFile(Y, resultFilePath, "Y");
            printVector(Y);

            writeResultsToFile(MA, resultFilePath, "MA");
            printMatrix(MA);
        }

        writeExecutionTimesToFile("results_lab2_wo_tsx/times.txt", times);

        System.out.println("Size,Times(ns)");
        for (Integer size : times.keySet()) {
            System.out.println(size + "," + times.get(size));
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

    public static void writeExecutionTimesToFile(String filePath, Map<Integer, Long> times) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            writer.write("Size,Times(ns)\n");
            for (Integer size : times.keySet()) {
                writer.write(size + "," + times.get(size) + "\n");
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

    public static void loadData(String filePath, int size) {
        MT = new double[size][size];
        MZ = new double[size][size];
        MA = new double[size][size];
        MT_plus_MZ = new double[size][size];
        MT_mlt_MT_plus_MZ = new double[size][size];
        MZ_mlt_MT = new double[size][size];
        MT_D = new double[size];
        B = new double[size];
        D = new double[size];
        Y = new double[size];

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

    static class RunnableThread implements Runnable {
        private final int start;
        private final int end;

        public RunnableThread(int start, int end) {
            this.start = start;
            this.end = end;
        }

        @Override
        public void run() {
            // Y = MT*D + max(B)*D
            double localMax = Double.NEGATIVE_INFINITY;
            for (int i = start; i <= end; i++) {
                localMax = Math.max(localMax, B[i]);
            }
            lockB_max.lock();
            try {
                b_max = Math.max(b_max, localMax);
            } finally {
                lockB_max.unlock();
            }
            
            // barrier to wait for b_max
            try {
                barrier.await();
            } catch (Exception e) {
                e.printStackTrace();
            }

            multiplyPartMatrixVectorKahan(MT, D, MT_D, start, end, lockMT_D);
            for (int i = start; i <= end; i++) {
                double result = MT_D[i] + (b_max * D[i]);
                lockY.lock();
                try {
                    Y[i] = result;
                } finally {
                    lockY.unlock();
                }
            }
            
            // MА = MT*(MT+MZ) - MZ*MT
            for (int i = 0; i < MT.length; i++) {
                for (int j = start; j <= end; j++) {
                    double result = MT[i][j] + MZ[i][j];
                    lockMT_plus_MZ.lock();
                    try {
                        MT_plus_MZ[i][j] = result;
                    } finally {
                        lockMT_plus_MZ.unlock();
                    }
                }
            }
            multiplyMatricesPartKahan(MT, MT_plus_MZ, MT_mlt_MT_plus_MZ, start, end, lockMT_mlt_MT_plus_MZ);
            multiplyMatricesPartKahan(MZ, MT, MZ_mlt_MT, start, end, lockMZ_mlt_MT);
            for (int i = 0; i < MT.length; i++) {
                for (int j = start; j <= end; j++) {
                    double result = MT_mlt_MT_plus_MZ[i][j] - MZ_mlt_MT[i][j];
                    lockMA.lock();
                    try {
                        MA[i][j] = result;
                    } finally {
                        lockMA.unlock();
                    }
                }
            }
        }

        private void multiplyPartMatrixVectorKahan(double[][] matrix, double[] vector, double[] result, int start, int end, Lock currentLock) {
            int size = vector.length;
    
            for (int i = start; i <= end; i++) {
                double sum = 0.0;
                double c = 0.0;
                for (int j = 0; j < size; j++) {
                    double y = matrix[i][j] * vector[j] - c;
                    double t = sum + y;
                    c = (t - sum) - y;
                    sum = t;
                }
    
                currentLock.lock();
                try {
                    result[i] = sum;
                } finally {
                    currentLock.unlock();
                }
            }
        }
    
        private void multiplyMatricesPartKahan(double[][] A, double[][] B, double[][] result, int start, int end, Lock currentLock) {
            for (int i = 0; i < B.length; i++) {
                for (int j = start; j <= end; j++) {
                    double sum = 0.0;
                    double c = 0.0;
                    for (int k = 0; k < B.length; k++) {
                        double y = A[i][k] * B[k][j] - c;
                        double t = sum + y;
                        c = (t - sum) - y;
                        sum = t;
                    }

                    currentLock.lock();
                    try {
                        result[i][j] = sum;
                    } finally {
                        currentLock.unlock();
                    }
                }
            }
        }
    }
}
