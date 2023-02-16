package Tetris.Graphics;

import Tetris.Operation.GameWindow;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import static javax.swing.LayoutStyle.ComponentPlacement.RELATED;

/**
 * Legmagasabb eredmények panel
 */
public class HighScores extends JPanel implements MouseListener {

    GameWindow gameWindow;
    JPanel scoresPanel;
    JLabel[] score;
    JLabel title;
    JButton back;
    Font newFont;

    /**
     * létrehozza a panelt és beolvassa a legmagasabb eredményeket
     * @param window az ablak amelyen a játék fut
     */
    public HighScores(GameWindow window){
    gameWindow = window;
    initialize();

    // Eredmények betöltése
    try {
        load();
    } catch (IOException e) {
        e.printStackTrace();
    }

    //Elrendezés
    GroupLayout layout = new GroupLayout(this);
    setLayout(layout);
    // Automatikus térköz
    layout.setAutoCreateGaps(true);
    // Vízszintes elrendezés
    layout.setHorizontalGroup(layout.createSequentialGroup()
            .addGroup(layout.createParallelGroup(GroupLayout.Alignment.CENTER)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(back,
                            GroupLayout.DEFAULT_SIZE, 70,
                            GroupLayout.PREFERRED_SIZE))
                        .addComponent(title,
                            GroupLayout.DEFAULT_SIZE, gameWindow.getWidth(),
                            GroupLayout.PREFERRED_SIZE)
                        .addComponent(scoresPanel))
    );
    // Függőleges elrendezés
    layout.setVerticalGroup(layout.createSequentialGroup()
            .addPreferredGap(RELATED, 5, 5)
            .addComponent(back,
                    GroupLayout.DEFAULT_SIZE, 70,
                    GroupLayout.PREFERRED_SIZE)
            .addPreferredGap(RELATED, 20, 20)
            .addComponent(title,
                    GroupLayout.DEFAULT_SIZE, 80,
                    GroupLayout.PREFERRED_SIZE)
            .addPreferredGap(RELATED, 40, 40)
            .addComponent(scoresPanel)
    );

    // Belső panel elrendezése
    GroupLayout layout2 = new GroupLayout(scoresPanel);
    scoresPanel.setLayout(layout2);
        // Vízszintes elrendezés
        layout2.setHorizontalGroup(layout2.createSequentialGroup()
                        .addGroup(layout2.createParallelGroup(GroupLayout.Alignment.CENTER)
                                .addComponent(score[0])
                                .addComponent(score[1])
                                .addComponent(score[2])
                                .addComponent(score[3])
                                .addComponent(score[4]))
        );
        // Függőleges elrendezés
        layout2.setVerticalGroup(layout2.createSequentialGroup()
                .addComponent(score[0])
                .addComponent(score[1])
                .addComponent(score[2])
                .addComponent(score[3])
                .addComponent(score[4])
        );

}

    /**
     * Inicializálja az adatokat
     */
    protected void initialize(){
        setBackground(Color.BLACK);

        title = new JLabel("High Scores");
        scoresPanel = new JPanel();
        try {
            newFont = getSpecialFont();
        }
        catch (IOException | FontFormatException exception){
            newFont = Font.getFont("Arial Bold");
        }

        title.setHorizontalAlignment(SwingConstants.CENTER);
        title.setBackground(Color.BLACK);
        title.setForeground(Color.RED);
        title.setBorder(null);
        title.setFont(newFont);

        back = new JButton();
        back.addMouseListener(this);
        try {
            Image img = ImageIO.read(new File("./res/backarrow.jpg"));
            back.setIcon(new ImageIcon(img));
            back.setHorizontalAlignment(SwingConstants.CENTER);
            back.setBackground(Color.BLACK);
            back.setBorder(null);
        } catch (Exception ex) {
            ex.printStackTrace();
        }


        scoresPanel.setBackground(Color.BLACK);
        score = new JLabel[5];
    }

    /**
     * Beolvassa az eredményeket
     * @throws IOException Nem sikerült megnyitni a fájlt
     */
    public void load() throws IOException{
        //beolvasás
        FileReader fr = new FileReader("./res/HighScores.txt");
        BufferedReader br = new BufferedReader(fr);
        int[] scores = new int[5];
        for (int i = 0; i < 5; i++){
            String line = br.readLine();
            if (line == null) break;
            scores[i] = Integer.parseInt(line);
        }
        br.close();
        // Eredmények betöltése
        for(int i = 0; i < 5; i++){
            score[i] = new JLabel(String.valueOf(scores[i]));
            score[i].setHorizontalAlignment(SwingConstants.CENTER);
            score[i].setBackground(Color.BLACK);
            score[i].setForeground(Color.BLUE);
            score[i].setBorder(null);
            score[i].setFont(newFont);
        }
    }

    /**
     * Betűtípus betöltése
     * @return betűtípus
     */
    public Font getSpecialFont() throws FontFormatException, IOException{
        Font font;
        // Fájlból olvasás és inicializálás
        font = Font.createFont(Font.TRUETYPE_FONT, new File("./res/Font.ttf")).deriveFont(40f);
        return font;
    }

    @Override
    public void mouseClicked(MouseEvent e) {
    }

    @Override
    public void mousePressed(MouseEvent e) {
    }

    @Override
    public void mouseReleased(MouseEvent e) {
            gameWindow.mainMenu(this);
    }

    @Override
    public void mouseEntered(MouseEvent e) {
        JButton button = (JButton) e.getComponent();
        button.setBackground(Color.WHITE);
    }

    @Override
    public void mouseExited(MouseEvent e) {
        JButton button = (JButton) e.getComponent();
        button.setBackground(Color.BLACK);
    }
}
