/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package dataanalysis;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JFrame;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.category.StackedBarRenderer;
import org.jfree.chart.renderer.category.StandardBarPainter;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.category.DefaultCategoryDataset;

/**
 *
 * @author yfa041
 */
public class BarChartView {

    public BarChartView(String title, double[][] data) {
        try {
            // Create dataset
            Thread t = new Thread(() -> {

                List<Double> ratios = calculateQuartileRatios(data);
                CategoryDataset dataset = createDataset(data);
                String updatedTitle = title + "-" + ratios;

                // Create chart
                JFreeChart chart = createChart(dataset, updatedTitle);

                // Customize the chart
                CategoryPlot plot = (CategoryPlot) chart.getPlot();
                StackedBarRenderer renderer = new StackedBarRenderer();
                chart.setBackgroundPaint(Color.WHITE);
                plot.setRenderer(renderer);

                // Customize renderer
                renderer.setSeriesPaint(1, Color.GRAY.brighter());
                renderer.setSeriesPaint(0, Color.GREEN.darker().darker());

                // Remove bar shadows
                renderer.setShadowVisible(false);
                plot.setBackgroundPaint(Color.WHITE);

                // Remove any reflections or glossy effects by making sure no effects are applied
                renderer.setBarPainter(new StandardBarPainter());
                // Create Panel
                ChartPanel panel = new ChartPanel(chart);
                JFrame frame = new JFrame(updatedTitle);
                frame.setSize(500, 500);
                frame.setContentPane(panel);
                frame.setVisible(true);
                frame.repaint();
            });
            t.start();
            t.join();
        } catch (InterruptedException ex) {
            Logger.getLogger(BarChartView.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    private List<Double> calculateQuartileRatios(double[][] data) {

        double totalCount = 0;
        List<Double> quartileRatios = new ArrayList<>();
        for (double[] qData : data) {
            totalCount += qData[0];
            quartileRatios.add(qData[0]);
        }
        double quat_3_sum = 0;
        for (int j = 0; j < quartileRatios.size() - 1; j++) {
            double ratio = quartileRatios.get(j) / totalCount;
            ratio = Math.round(ratio * 100.0) / 100.0;
            quartileRatios.set(j, ratio);
            quat_3_sum += ratio;
        }
        double ratio = 1.0 - quat_3_sum;
        ratio = Math.round(ratio * 100.0) / 100.0;
        quartileRatios.set(quartileRatios.size() - 1, ratio);
        
        quartileRatios.add(totalCount);
        return quartileRatios;

    }

    private CategoryDataset createDataset(double[][] data) {
        DefaultCategoryDataset dataset = new DefaultCategoryDataset();
        int i = 1;
        for (double[] qData : data) {
            dataset.addValue(qData[0], "Identified", "Q" + i);
            dataset.addValue(qData[1], "UnIdentified", "Q" + i);
            i++;
        }
        return dataset;
    }

    private JFreeChart createChart(CategoryDataset dataset, String title) {
        return ChartFactory.createStackedBarChart(
                title, // Chart title
                "Category", // X-Axis Label
                "Value", // Y-Axis Label
                dataset, // Dataset
                PlotOrientation.VERTICAL, // Orientation
                true, // Include legend
                true, // Tooltips
                false // URLs?
        );
    }
}
