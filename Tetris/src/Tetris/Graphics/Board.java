package Tetris.Graphics;

import javax.swing.*;
import java.awt.*;

/**
 * Megjeleníti a táblát
 */
public class Board extends JPanel {

    private static final int WIDTH = 26 * 14;
    private static final int HEIGHT = 26 * 25;
    private Color[][] well;
    JLabel pauseText = new JLabel("PAUSED");
    private boolean paused;

    /**
     * Létrehoz egy tábla panelt
     */
    public Board() {
        well = new Color[14][24];
        for (int i = 0; i < 14; i++) {
            for (int j = 0; j < 24; j++) {
                if (i == 0 || i == 1 || i == 12 || i == 13 || j == 22 || j == 23) {
                    well[i][j] = Color.GRAY;
                } else {
                    well[i][j] = Color.BLACK;
                }
            }
        }
    }

    /**
     * Kirajzolja a táblát
     * @param g Graphics
     */
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);

        g.fillRect(0, 0, WIDTH, HEIGHT);
        for (int i = 0; i < 14; i++) {
            for (int j = 0; j < 24; j++) {
                g.setColor(well[i][j]);
                g.fillRect(26 * i, 26 * j, 25, 25);
            }
        }

        if (paused) {
            add(pauseText);
            pauseText.setForeground(Color.WHITE);
            pauseText.setFont(pauseText.getFont().deriveFont(50f));
            validate();
        } else remove(pauseText);

    }

    /**
     * Újrarajzolja a táblát
     * @param w tábla
     */
    public void paintWell(Color[][] w) {
        well = w;
        repaint();
    }

    /**
     * Aktiválja/Deaktiválja a szünet panelt
     * @param status panel állapota
     */
    public void pausePanel(boolean status) {
        paused = status;
    }
}
