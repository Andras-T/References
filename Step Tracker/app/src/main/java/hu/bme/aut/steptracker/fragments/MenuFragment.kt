package hu.bme.aut.steptracker.fragments

import android.os.Bundle
import androidx.fragment.app.Fragment
import android.view.LayoutInflater
import android.view.View
import android.view.ViewGroup
import androidx.navigation.fragment.findNavController
import hu.bme.aut.steptracker.R
import hu.bme.aut.steptracker.databinding.FragmentMenuBinding

class MenuFragment : Fragment() {
    lateinit var binding : FragmentMenuBinding

    override fun onCreateView(inflater: LayoutInflater, container: ViewGroup?, savedInstanceState: Bundle?): View {
        super.onCreateView(inflater, container, savedInstanceState)
        binding = FragmentMenuBinding.inflate(inflater, container, false)
        return binding.root
    }

    override fun onViewCreated(view: View, savedInstanceState: Bundle?) {
        super.onViewCreated(view, savedInstanceState)

        binding.trackerBtn.setOnClickListener {
            findNavController().navigate(R.id.action_menuFragment_to_infoFragment)
        }
        binding.dataBtn.setOnClickListener {
            findNavController().navigate(R.id.action_menuFragment_to_dataFragment)
        }
        binding.achievementsBtn.setOnClickListener {
            findNavController().navigate(R.id.action_menuFragment_to_achievementFragment)
        }
        binding.settingsBtn.setOnClickListener {
            findNavController().navigate(R.id.action_menuFragment_to_settingsFragment)
        }
    }
}