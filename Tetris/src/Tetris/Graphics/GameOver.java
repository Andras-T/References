package Tetris.Graphics;

import Tetris.Operation.GameWindow;
import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import java.io.IOException;
import static javax.swing.LayoutStyle.ComponentPlacement.RELATED;

/**
 * Játék vége panel
 */
public class GameOver extends JPanel implements MouseListener {

    GameWindow gameWindow;
    JLabel over, scoreLabel;
    JButton back;

    /**
     * Panel létrehozása és komponensek beállítása
     * @param window az ablak amin megjelenik
     * @param score elért eredmény
     */
    public GameOver(GameWindow window, int score){
        gameWindow = window;
        initialize(score);

        // Elrendezés
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
                        .addGroup(layout.createSequentialGroup()
                            .addComponent(over,
                                GroupLayout.DEFAULT_SIZE, gameWindow.getWidth(),
                                GroupLayout.PREFERRED_SIZE))
                       .addComponent(scoreLabel,
                                GroupLayout.DEFAULT_SIZE, gameWindow.getWidth(),
                                GroupLayout.PREFERRED_SIZE))
        );
        // Függőleges elrendezés
        layout.setVerticalGroup(layout.createSequentialGroup()
                .addPreferredGap(RELATED, 5, 5)
                .addComponent(back,
                        GroupLayout.DEFAULT_SIZE, 70,
                        GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(RELATED, 20, 20)
                .addComponent(over,
                        GroupLayout.DEFAULT_SIZE, 80,
                        GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(RELATED, 20, 20)
                .addComponent(scoreLabel,
                        GroupLayout.DEFAULT_SIZE, 80,
                        GroupLayout.PREFERRED_SIZE)
        );
    }

    /**
     * Inicializálja az adatokat
     * @param score elért eredmény
     */
    protected void initialize(int score){
        setBackground(Color.BLACK);

        over = new JLabel("GAME OVER!");
        over.setForeground(Color.RED);
        over.setHorizontalAlignment(SwingConstants.CENTER);

        scoreLabel = new JLabel("Points: " + score);
        scoreLabel.setForeground(Color.WHITE);
        scoreLabel.setHorizontalAlignment(SwingConstants.CENTER);

        back = new JButton("Return to the menu");
        back.setHorizontalAlignment(SwingConstants.CENTER);
        back.setBackground(Color.BLACK);
        back.setForeground(Color.GREEN);
        back.setBorder(null);
        back.addMouseListener(this);

        // Betűtípus beállítása
        Font newFont = Font.getFont("Arial Bold");
        try {
            newFont = getSpecialFont();
            over.setFont(newFont.deriveFont(50f));
            scoreLabel.setFont(newFont.deriveFont(35f));
            back.setFont(newFont.deriveFont(25f));
        } catch (IOException | FontFormatException e) {
            over.setFont(over.getFont().deriveFont(50f));
            scoreLabel.setFont(scoreLabel.getFont().deriveFont(35f));
            back.setFont(back.getFont().deriveFont(25f));
        }
    }

    /**
     * Beolvas egy betűtípust egy fájlból
     * @return Betűtípus
     * @throws IOException Nem sikerült beolvasni a fájlt
     * @throws FontFormatException Hibás fájl tartalom
     */
    public Font getSpecialFont() throws IOException, FontFormatException{
        Font font;
        // Fájlból olvasás
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
