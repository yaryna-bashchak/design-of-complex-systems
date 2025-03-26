package lab3.src.main.java.lab3;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.Statement;
import java.sql.SQLException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import static org.multiverse.api.StmUtils.atomic;

public class LoadTest {
    public static void main(String[] args) {
        String url = "jdbc:mysql://localhost:3306/test_db";
        String user = "root";
        String password = "TempYabas7/";
        String outputFile = "htm_results.txt";
        int numThreads = 4;

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
            writer.write("Number of records,Time(ms)\n");

            for (int numRecordsPerThread = 100; numRecordsPerThread <= 1000; numRecordsPerThread += 100) {

                ExecutorService executor = Executors.newFixedThreadPool(numThreads);

                long startTime = System.currentTimeMillis();

                for (int i = 0; i < numThreads; i++) {
                    executor.submit(new InsertTask(url, user, password, numRecordsPerThread, i + 1));
                }

                executor.shutdown();
                while (!executor.isTerminated()) {
                }

                long endTime = System.currentTimeMillis();
                long duration = endTime - startTime;

                writer.write(numThreads * numRecordsPerThread + ", " + duration + "\n");

                System.out.println("Inserted " + (numThreads * numRecordsPerThread) + " records in " + duration
                        + " milliseconds.");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    static class InsertTask implements Runnable {
        private String url;
        private String user;
        private String password;
        private int numRecords;
        private int threadId;

        public InsertTask(String url, String user, String password, int numRecords, int threadId) {
            this.url = url;
            this.user = user;
            this.password = password;
            this.numRecords = numRecords;
            this.threadId = threadId;
        }

        @Override
        public void run() {
            try (Connection conn = DriverManager.getConnection(url, user, password)) {
                if (conn != null) {
                    Statement stmt = conn.createStatement();
                    String insertSqlTemplate = "INSERT INTO users(name, age) VALUES('Alice_%d_%d', 30);";

                    for (int i = 0; i < numRecords; i++) {
                        String insertSql = String.format(insertSqlTemplate, threadId, i);
                        
                        // comment one of the following blocks
                        ///* HTM
                        stmt.executeUpdate(insertSql);
                        //*/

                        ///* STM
                        atomic(() -> {
                            try {
                                stmt.executeUpdate(insertSql);
                            } catch (SQLException e) {
                                e.printStackTrace();
                            }
                        });
                        //*/
                    }
                }
            } catch (SQLException e) {
                System.out.println("Error: " + e.getMessage());
            }
        }
    }
}
