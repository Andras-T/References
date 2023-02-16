package Tetris.Graphics;

import Tetris.Operation.GameWindow;
import javax.swing.*;
import javax.swing.border.LineBorder;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.font.TextAttribute;
import java.io.File;
import java.io.IOException;
import java.util.Map;
import static javax.swing.LayoutStyle.ComponentPlacement.RELATED;

/**
 * Főmenü panel
 */
public class MainMenu extends JPanel implements MouseListener {

    JButton play, scores, quit;
    JTextField gameName, creator;
    private Font buttonFont;
    private final Font titleFont;
    private final LineBorder border;
    private static final int SPACING = 20;
    private final GameWindow gameWindow;

    final String buttonName1 = "Play";
    final String buttonName2 = "High Scores";
    final String buttonName3 = "Quit";

    /**
     * Főmenü panel létrehozása
     * @param window az ablak amelyen a játék fut
     */
    public MainMenu(GameWindow window){
        gameWindow = window;
        setBackground(Color.BLACK);

        final String Titel = "TETRIS";
        border = new LineBorder(Color.BLACK,1);
        play = new JButton(buttonName1);
        scores = new JButton(buttonName2);
        quit = new JButton(buttonName3);
        gameName = new JTextField(Titel);
        creator = new JTextField("Tóth András");
        GroupLayout layout = new GroupLayout(this);

        play.setName(buttonName1);
        scores.setName(buttonName2);
        quit.setName(buttonName3);

        play.addMouseListener(this);
        scores.addMouseListener(this);
        quit.addMouseListener(this);

        // Betűtípusok beállítása
        try {
            buttonFont = getSpecialFont();
        } catch (IOException | FontFormatException e) {
            buttonFont = Font.getFont("Arial Bold");
        }
        titleFont = buttonFont.deriveFont(80f);

        // Elrendezés beállítása
        setLayout(layout);
        // Automatikus térköz
        layout.setAutoCreateGaps(true);
        // Vízszintes elrendezés
        layout.setHorizontalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.CENTER)
                        .addComponent(gameName,
                        GroupLayout.DEFAULT_SIZE, gameWindow.getWidth(),
                        GroupLayout.PREFERRED_SIZE)
                        .addComponent(play)
                        .addComponent(scores)
                        .addComponent(quit)
                        .addComponent(creator, GroupLayout.DEFAULT_SIZE, 150,
                                GroupLayout.PREFERRED_SIZE))
                );
        // Függőleges elrendezés
        layout.setVerticalGroup(layout.createSequentialGroup()
                .addComponent(gameName,
                        GroupLayout.DEFAULT_SIZE, 50,
                        GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(RELATED, 70, 70)
                .addComponent(play)
                .addPreferredGap(RELATED, SPACING, SPACING)
                .addComponent(scores)
                .addPreferredGap(RELATED, SPACING, SPACING)
                .addComponent(quit)
                .addPreferredGap(RELATED,
                        GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(creator, GroupLayout.DEFAULT_SIZE, 50,
                        GroupLayout.PREFERRED_SIZE)
);

    }

    /**
     * Fájlból beolvas egy betűtípust
     */
    public Font getSpecialFont() throws IOException, FontFormatException {
        Font font;
        //fájlból olvasás
        font = Font.createFont(Font.TRUETYPE_FONT, new File("./res/Font.ttf")).deriveFont(40f);
        return font;
    }

    /**
     * Kirajzolja  a komponenseket
     * @param g Graphics
     */
    protected void paintComponent(Graphics g){
        super.paintComponent(g);
        g.setColor(Color.BLACK);
        g.fillRect(0, 0, 26*12, 26*23);
        paintButton(play, Color.ORANGE);
        paintButton(scores, Color.CYAN);
        paintButton(quit, Color.GREEN);

        gameName.setEditable(false);
        gameName.setHorizontalAlignment(SwingConstants.CENTER);
        gameName.setBackground(Color.BLACK);
        gameName.setForeground(Color.RED);
        gameName.setBorder(border);
        // Aláhúzás
        Map attributes = titleFont.getAttributes();
        attributes.put(TextAttribute.UNDERLINE, TextAttribute.UNDERLINE_ON);
        gameName.setFont(titleFont.deriveFont(attributes));

        creator.setEditable(false);
        creator.setHorizontalAlignment(SwingConstants.CENTER);
        creator.setBackground(Color.BLACK);
        creator.setBorder(border);
    }

    /**
     * Beállítja a gombok kinézetét
     * @param button Az adott gomb
     * @param color A gomb színe
     */
    protected void paintButton(JButton button, Color color){
        button.setForeground(color);
        button.setBackground(Color.BLACK);
        button.setFont(buttonFont);
        button.setBorder(border);
    }

    /*--------------------------Irányítás----------------------------------------*/

    @Override
    public void mouseClicked(MouseEvent e) {}

    @Override
    public void mousePressed(MouseEvent e) {}

    @Override
    public void mouseReleased(MouseEvent e) {
        if(buttonName1.equals(e.getComponent().getName()))
            gameWindow.start();
        else if(buttonName2.equals(e.getComponent().getName()))
            gameWindow.highScores();
        else if(buttonName3.equals(e.getComponent().getName())){
            gameWindow.close();
        }
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